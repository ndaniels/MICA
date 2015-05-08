package main

import (
	"math"
	"runtime"
	"sync"
	//	"time"

	"github.com/ndaniels/neutronium"
)

type alignPool struct {
	db            *neutronium.DB
	jobs          chan *alignJob
	results       chan *seqComparison
	bestMatch     *seqComparison
	jobWG         *sync.WaitGroup
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
func startCompressWorkers(db *neutronium.DB, seedTable *neutronium.SeedTable, mems []*memory) alignPool {
	jobWG := &sync.WaitGroup{}
	jobs := make(chan *alignJob, 200)
	results := make(chan *seqComparison, 10*1000)

	pool := alignPool{
		db:            db,
		jobs:          jobs,
		results:       results,
		bestMatch:     nil,
		jobWG:         jobWG,
		closed:        false,
		allJobsLoaded: false,
		seedTable:     seedTable,
	}
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
		jobWG.Add(1)
		go pool.aligner(mems[i])
	}

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

func (pool *alignPool) aligner(mem *memory) {
	maxRadius := pool.db.DBConf.MaxClusterRadius
	for job := range pool.jobs {
		comp := compareSeqs(maxRadius, job.corSeqId, job.orgSeqId, job.corSeq, job.orgSeq, pool.seedTable, mem)
		if comp.distance <= maxRadius {
			pool.results <- comp
		}
	}
	pool.jobWG.Done()
}

func (pool *alignPool) emptyResultsQueue() {
	for comp := range pool.results {
		if pool.bestMatch != nil {
			if math.Abs(comp.distance-pool.bestMatch.distance) < 0.000001 &&
				len(comp.corSeq.Residues) > len(pool.bestMatch.corSeq.Residues) {
				pool.bestMatch = comp
			} else if comp.distance < pool.bestMatch.distance {
				pool.bestMatch = comp
			}
		} else {
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

	// The bestMatch in pool should now actually be the best match, if it exists
	pool.emptyResultsQueue()
	bestMatch := pool.bestMatch

	if bestMatch == nil || bestMatch.distance > pool.db.DBConf.MaxClusterRadius {
		newComSeq := neutronium.NewCompressedSeq(pool.seqId, pool.seq.Name)
		addWithoutMatch(&newComSeq, pool.db.CoarseDB, pool.seqId, pool.seq, pool.seedTable)
	} else {
		addWithMatch(bestMatch.orgSeq, bestMatch.corSeq, bestMatch.alignment, bestMatch.orgSeqId, bestMatch.corSeqId)
	}

}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
func addWithoutMatch(cseq *neutronium.CompressedSeq, coarsedb *neutronium.CoarseDB, orgSeqId int, orgSub *neutronium.OriginalSeq, seedTable *neutronium.SeedTable) int {
	// Explicitly copy residues to avoid pinning memory.
	subCpy := make([]byte, len(orgSub.Residues))
	copy(subCpy, orgSub.Residues)

	corSeqId, corSeq := coarsedb.Add(subCpy)
	corSeq.AddLink(
		neutronium.NewLinkToCompressed(uint32(orgSeqId), 0, uint16(len(subCpy))))

	cseq.Add(
		neutronium.NewLinkToCoarseNoDiff(uint(corSeqId), 0, uint(len(subCpy))))

	seedTable.Add(corSeqId, corSeq)

	return corSeqId
}

func addWithMatch(oSeq *neutronium.OriginalSeq, corSeq *neutronium.CoarseSeq, alignment [2][]byte, oSeqId, corSeqId int) {

	comSeq := neutronium.NewCompressedSeq(oSeqId, oSeq.Name)

	corLen := corSeq.Len()

	coarseLink := neutronium.NewLinkToCoarse(uint(corSeqId), 0, uint(corLen), alignment)
	comSeq.Add(coarseLink)

	compressedLink := neutronium.NewLinkToCompressed(uint32(oSeqId), 0, uint16(corLen))
	corSeq.AddLink(compressedLink)
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
	minPossibleDistance := 1.0 - float64(min(cLen, oLen))/float64(max(cLen, oLen))

	if minPossibleDistance > matchThreshold {
		return &seqComparison{
			distance: 1,
			corSeqId: corSeqId,
			corSeq:   corSeq,
			orgSeqId: orgSeqId,
			orgSeq:   orgSeq,
		}
	}

	m := max(cLen, oLen)
	k := seedTable.SeedSize
	gaps := m - min(cLen, oLen)
	allowableMismatches := int((float64(m) * matchThreshold)) - gaps
	maxMissingMers := allowableMismatches * k
	numberMersWithoutMisses := m - k + 1
	matchingKmerLowerBound := numberMersWithoutMisses - maxMissingMers

	matchingKmers := 0
	for i := 0; i < orgSeq.Len()-k; i++ {
		kmer := orgSeq.Residues[i : i+k]
		if seedTable.Lookup(kmer, corSeqId) {
			matchingKmers++
			if matchingKmers >= matchingKmerLowerBound {
				break
			}
		} else if orgSeq.Len()-k-i < matchingKmerLowerBound {
			break
		}
	}

	if matchingKmers < matchingKmerLowerBound {
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
