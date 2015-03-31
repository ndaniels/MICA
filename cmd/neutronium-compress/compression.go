package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/ndaniels/neutronium"
)

type alignPool struct {
	db            *neutronium.DB
	jobs          chan *alignJob
	results       chan *seqComparison
	bestMatch     *seqComparison
	jobWG         *sync.WaitGroup
	recWG         *sync.WaitGroup
	closed        bool
	allJobsLoaded bool
	seqId         int
	seq           *neutronium.OriginalSeq
	seedTable     *neutronium.SeedTable
}

type alignJob struct {
	orgSeqId int
	orgSeq   *neutronium.OriginalSeq
	corSeqId int
	corSeq   *neutronium.CoarseSeq
}

// startCompressWorkers initializes a pool of compression workers.
//
// The compressPool returned can be used to compress sequences concurrently.
func startCompressWorkers(db *neutronium.DB, seedTable *neutronium.SeedTable) alignPool {
	jobWG := &sync.WaitGroup{}
	recWG := &sync.WaitGroup{}
	jobs := make(chan *alignJob, 200)
	results := make(chan *seqComparison, 200)
	pool := alignPool{
		db:            db,
		jobs:          jobs,
		results:       results,
		bestMatch:     nil,
		jobWG:         jobWG,
		recWG:         recWG,
		closed:        false,
		allJobsLoaded: false,
		seedTable:     seedTable,
	}
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)-12); i++ {
		jobWG.Add(1)
		go pool.aligner()
	}
	recWG.Add(1)
	go pool.receiver(pool.db.DBConf.MaxClusterRadius)
	return pool
}

func (pool *alignPool) align(id int, seq *neutronium.OriginalSeq) {

	pool.seqId = id
	pool.seq = seq

	coarseDB := pool.db.CoarseDB
	for cSeqId, cSeq := range coarseDB.Seqs {
		pool.jobs <- &alignJob{
			orgSeqId: id,
			orgSeq:   seq,
			corSeqId: cSeqId,
			corSeq:   cSeq,
		}
	}
	pool.allJobsLoaded = true
	close(pool.jobs)
}

func (pool *alignPool) aligner() {
	mem := newMemory()
	for job := range pool.jobs {
		comp := compareSeqs(pool.db.DBConf.MaxClusterRadius, job.corSeqId, job.orgSeqId, job.corSeq, job.orgSeq, pool.seedTable, mem)
		pool.results <- comp
	}
	pool.jobWG.Done()
}

func (pool *alignPool) receiver(maxRadius float64) {
	pool.emptyResultsQueue(maxRadius)
	pool.recWG.Done()
}

func (pool *alignPool) emptyResultsQueue(maxRadius float64) {
	for comp := range pool.results {
		if pool.bestMatch != nil {
			if math.Abs(comp.distance-pool.bestMatch.distance) < 0.000001 &&
				len(comp.corSeq.Residues) > len(pool.bestMatch.corSeq.Residues) {
				pool.bestMatch = comp
			} else if comp.distance < pool.bestMatch.distance &&
				comp.distance <= maxRadius {
				pool.bestMatch = comp
			}
		} else if comp.distance <= maxRadius {
			pool.bestMatch = comp
		}
	}
}

func (pool *alignPool) finishAndHandle() {
	if pool.closed {
		return
	}
	pool.closed = true
	pool.jobWG.Wait()
	close(pool.results)
	pool.recWG.Wait()
	// The bestMatch in pool should now actually be the best match, if it exists
	bestMatch := pool.bestMatch
	newComSeq := neutronium.NewCompressedSeq(pool.seqId, pool.seq.Name)
	if pool.bestMatch == nil {
		// Add the input sequence as a new representative
		addWithoutMatch(&newComSeq, pool.db.CoarseDB, pool.seqId, pool.seq, pool.seedTable)
	} else {
		// Link the input sequence to the match
		corLen := uint(len(bestMatch.corSeq.Residues))
		newComSeq.Add(neutronium.NewLinkToCoarse(
			uint(bestMatch.corSeqId), 0, corLen, bestMatch.alignment))
		bestMatch.corSeq.AddLink(neutronium.NewLinkToCompressed(
			uint32(bestMatch.orgSeqId), 0, uint16(corLen)))
	}

}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
func addWithoutMatch(cseq *neutronium.CompressedSeq, coarsedb *neutronium.CoarseDB, orgSeqId int, orgSub *neutronium.OriginalSeq, seedTable *neutronium.SeedTable) {

	// Explicitly copy residues to avoid pinning memory.
	subCpy := make([]byte, len(orgSub.Residues))
	copy(subCpy, orgSub.Residues)

	corSeqId, corSeq := coarsedb.Add(subCpy)
	corSeq.AddLink(
		neutronium.NewLinkToCompressed(uint32(orgSeqId), 0, uint16(len(subCpy))))

	cseq.Add(
		neutronium.NewLinkToCoarseNoDiff(uint(corSeqId), 0, uint(len(subCpy))))

	seedTable.Add(corSeqId, corSeq)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

type seqComparison struct {
	distance  float64
	alignment [2][]byte
	corSeqId  int
	corSeq    *neutronium.CoarseSeq
	orgSeqId  int
	orgSeq    *neutronium.OriginalSeq
}

func compareSeqs(matchThreshold float64, corSeqId, orgSeqId int, corSeq *neutronium.CoarseSeq, orgSeq *neutronium.OriginalSeq, seedTable *neutronium.SeedTable, mem *memory) *seqComparison {
	cLen := len(corSeq.Residues)
	oLen := len(orgSeq.Residues)
	minDistance := 1.0 - float64(min(cLen, oLen))/float64(max(cLen, oLen))

	if minDistance > matchThreshold {
		return &seqComparison{
			distance: 1,
			corSeqId: corSeqId,
			corSeq:   corSeq,
			orgSeqId: orgSeqId,
			orgSeq:   orgSeq,
		}
	}

	matchingKmers := 0
	for i := 0; i < orgSeq.Len()-seedTable.SeedSize; i++ {
		kmer := orgSeq.Residues[i : i+seedTable.SeedSize]
		if seedTable.Lookup(kmer, corSeqId) {
			matchingKmers++
		}
	}

	if matchingKmers < 1000 {
		return &seqComparison{
			distance: 1,
			corSeqId: corSeqId,
			corSeq:   corSeq,
			orgSeqId: orgSeqId,
			orgSeq:   orgSeq,
		}
	}

	alignment := nwAlign(corSeq.Residues, orgSeq.Residues, mem)
	distance := alignmentDistance(alignment)

	return &seqComparison{
		distance:  distance,
		alignment: alignment,
		corSeqId:  corSeqId,
		corSeq:    corSeq,
		orgSeqId:  orgSeqId,
		orgSeq:    orgSeq,
	}

}
