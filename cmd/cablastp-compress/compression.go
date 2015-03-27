package main

import (
	"bytes"
	"runtime"
	"sync"

	"github.com/BurntSushi/cablastp"
)

type alignPool struct {
	db     *cablastp.DB
	jobs   chan alignJob
	wg     *sync.WaitGroup
	closed bool
}

type alignJob struct {
	orgSeqId int
	orgSeq   *cablastp.OriginalSeq
	corSeqId int
	corSeq   *cablastp.CoarseSeq
}

// startCompressWorkers initializes a pool of compression workers.
//
// The compressPool returned can be used to compress sequences concurrently.
func startCompressWorkers(db *cablastp.DB) alignPool {
	wg := &sync.WaitGroup{}
	jobs := make(chan alignJob, 200)
	pool := alignPool{
		db:     db,
		jobs:   jobs,
		wg:     wg,
		closed: false,
	}
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
		wg.Add(1)
		go pool.worker()
	}
	return pool
}

func (pool compressPool) align(id int, seq *cablastp.OriginalSeq) {
	// Still needs seed table optimization
	for corSeqId = 0; corSeqId < db.CoarseDB.NumSequences(); corSeqId++ {
		pool.jobs <- alignJob{
			orgSeqId: id,
			orgSeq:   seq,
			corSeqId: corSeqId,
			corSeq: db.CoarseDb.CoarseSeqGet(corSeqId)
		}
	}
}


func (pool compressPool) worker() {
	mem := newMemory()
	for job := range pool.jobs {
		comp := compareSeqs( pool.db.DBConf.MaxClusterRadius, job.corSeq, job.orgSeq, mem)
		// What do I do with the output from compareSeqs?
	}
	pool.wg.Done()
}

func (pool *compressPool) done() {
	if pool.closed {
		return
	}
	pool.closed = true
	close(pool.jobs)
	pool.wg.Wait()
}


// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
func addWithoutMatch(cseq *cablastp.CompressedSeq,
	coarsedb *cablastp.CoarseDB, orgSeqId int, orgSub *cablastp.OriginalSeq) {

	// Explicitly copy residues to avoid pinning memory.
	subCpy := make([]byte, len(orgSub.Residues))
	copy(subCpy, orgSub.Residues)

	corSeqId, corSeq := coarsedb.Add(subCpy)
	corSeq.AddLink(
		cablastp.NewLinkToCompressed(uint32(orgSeqId), 0, uint16(len(subCpy))))

	cseq.Add(
		cablastp.NewLinkToCoarseNoDiff(uint(corSeqId), 0, uint(len(subCpy))))
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
	corSeq *cablastp.CoarseSeq
	orgSeq *cablastp.OriginalSeq
}

func compareSeqs(matchThreshold float64, cseq *cablastp.CoarseSeq, oseq *cablastp.OriginalSeq, mem *memory) seqComparison {
	cLen := len(cseq.Residues)
	oLen := len(oseq.Residues)
	minDistance := 1 - min(cLen, oLen)/max(cLen, oLen)

	if minDistance > matchThreshold {
		return seqComparison{
			distance:  1,
			alignment: nil,
			corSeq: cseq,
			orgSeq: oseq,
		}
	}

	alignment := nwAlign(cseq.Residues, oseq.Residues, mem)
	distance := alignmentDistance(alignment)

	return seqComparison{
		distance:  distance,
		alignment: alignment,
		corSeq: cseq,
		orgSeq: oseq,
	}

}
