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
	wg            *sync.WaitGroup
	closed        bool
	allJobsLoaded bool
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
func startCompressWorkers(db *neutronium.DB) alignPool {
	wg := &sync.WaitGroup{}
	jobs := make(chan *alignJob, 200)
	results := make(chan *seqComparison, 200)
	pool := alignPool{
		db:            db,
		jobs:          jobs,
		results:       results,
		bestMatch:     nil,
		wg:            wg,
		closed:        false,
		allJobsLoaded: false,
	}
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
		wg.Add(1)
		go pool.aligner()
	}
	wg.Add(1)
	go pool.receiver(pool.db.DBConf.MaxClusterRadius)
	return pool
}

func (pool *alignPool) align(id int, seq *neutronium.OriginalSeq) {
	// Still needs seed table optimization
	coarseDB := pool.db.CoarseDB
	for corSeqId := 0; corSeqId < coarseDB.NumSequences(); corSeqId++ {
		pool.jobs <- &alignJob{
			orgSeqId: id,
			orgSeq:   seq,
			corSeqId: corSeqId,
			corSeq:   coarseDB.CoarseSeqGet(uint(corSeqId)),
		}
	}
	pool.allJobsLoaded = true
}

func (pool *alignPool) aligner() {
	mem := newMemory()
	for job := range pool.jobs {
		comp := compareSeqs(pool.db.DBConf.MaxClusterRadius, job.corSeqId, job.orgSeqId, job.corSeq, job.orgSeq, mem)
		pool.results <- &comp
	}
	pool.wg.Done()
}

func (pool *alignPool) receiver(maxRadius float64) {
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
	if !pool.allJobsLoaded {
		pool.wg.Add(1)
		pool.receiver(maxRadius)
	}
	pool.wg.Done()
}

func (pool *alignPool) finishAndHandle() {
	if pool.closed {
		return
	}
	pool.closed = true
	close(pool.jobs)
	pool.wg.Wait()
	pool.wg.Add(1) // Avoid a panic when the receiver decrements the wait group
	pool.receiver(pool.db.DBConf.MaxClusterRadius)
	close(pool.results)
	// The bestMatch in pool should now actually be the best match
	newComSeq := neutronium.NewCompressedSeq(pool.bestMatch.orgSeqId, pool.bestMatch.orgSeq.Name)
	if pool.bestMatch == nil {
		addWithoutMatch(&newComSeq, pool.db.CoarseDB, pool.bestMatch.orgSeqId, pool.bestMatch.orgSeq)
	} else {
		corLen := uint(len(pool.bestMatch.corSeq.Residues))
		newComSeq.Add(neutronium.NewLinkToCoarse(
			uint(pool.bestMatch.corSeqId), 0, corLen, pool.bestMatch.alignment))
		pool.bestMatch.corSeq.AddLink(neutronium.NewLinkToCompressed(
			uint32(pool.bestMatch.orgSeqId), 0, uint16(corLen)))
	}

}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
func addWithoutMatch(cseq *neutronium.CompressedSeq,
	coarsedb *neutronium.CoarseDB, orgSeqId int, orgSub *neutronium.OriginalSeq) {

	// Explicitly copy residues to avoid pinning memory.
	subCpy := make([]byte, len(orgSub.Residues))
	copy(subCpy, orgSub.Residues)

	corSeqId, corSeq := coarsedb.Add(subCpy)
	corSeq.AddLink(
		neutronium.NewLinkToCompressed(uint32(orgSeqId), 0, uint16(len(subCpy))))

	cseq.Add(
		neutronium.NewLinkToCoarseNoDiff(uint(corSeqId), 0, uint(len(subCpy))))
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

func compareSeqs(matchThreshold float64, corSeqId, orgSeqId int, corSeq *neutronium.CoarseSeq, orgSeq *neutronium.OriginalSeq, mem *memory) seqComparison {
	cLen := len(corSeq.Residues)
	oLen := len(orgSeq.Residues)
	minDistance := 1.0 - float64(min(cLen, oLen))/float64(max(cLen, oLen))

	if minDistance > matchThreshold {
		return seqComparison{
			distance: 1,
			corSeqId: corSeqId,
			corSeq:   corSeq,
			orgSeqId: orgSeqId,
			orgSeq:   orgSeq,
		}
	}

	alignment := nwAlign(corSeq.Residues, orgSeq.Residues, mem)
	distance := alignmentDistance(alignment)

	return seqComparison{
		distance:  distance,
		alignment: alignment,
		corSeqId:  corSeqId,
		corSeq:    corSeq,
		orgSeqId:  orgSeqId,
		orgSeq:    orgSeq,
	}

}
