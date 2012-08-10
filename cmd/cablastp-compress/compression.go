package main

import (
	"bytes"
	"runtime"
	"sync"

	"github.com/BurntSushi/cablastp"
)

type compressPool struct {
	db   *cablastp.DB
	jobs chan compressJob
	wg   *sync.WaitGroup
}

type compressJob struct {
	orgSeqId int
	orgSeq   *cablastp.OriginalSeq
}

func startCompressWorkers(db *cablastp.DB) compressPool {
	wg := &sync.WaitGroup{}
	jobs := make(chan compressJob, 200)
	pool := compressPool{
		db:   db,
		jobs: jobs,
		wg:   wg,
	}
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)); i++ {
		wg.Add(1)
		go pool.worker()
	}
	return pool
}

func (pool compressPool) compress(id int, seq *cablastp.OriginalSeq) int {
	pool.jobs <- compressJob{
		orgSeqId: id,
		orgSeq:   seq,
	}
	return id + 1
}

func (pool compressPool) worker() {
	mem := newMemory()
	for job := range pool.jobs {
		comSeq := compress(pool.db, job.orgSeqId, job.orgSeq, mem)
		pool.db.ComDB.Write(comSeq)
	}
	pool.wg.Done()
}

func (pool compressPool) done() {
	close(pool.jobs)
	pool.wg.Wait()
}

func compress(db *cablastp.DB, orgSeqId int,
	orgSeq *cablastp.OriginalSeq, mem *memory) cablastp.CompressedSeq {

	var cseqExt, oseqExt []byte

	coarsedb := db.CoarseDB
	cseq := cablastp.NewCompressedSeq(orgSeqId, orgSeq.Name)
	mapSeedSize := coarsedb.Seeds.SeedSize
	extSeedSize := db.ExtSeedSize

	// Keep track of two pointers. 'current' refers to the residue index in the
	// original sequence that extension is currently originating from.
	// 'lastMatch' refers to the residue index of the *end* of the last match
	// with a reference sequence in the compressed database.
	lastMatch, current := 0, 0

	// Iterate through the original sequence a 'kmer' at a time.
	for current = 0; current < orgSeq.Len()-mapSeedSize-extSeedSize; current++ {
		kmer := orgSeq.Residues[current : current+mapSeedSize]
		seeds := coarsedb.Seeds.Lookup(kmer, &mem.seeds)
		if seeds == nil {
			continue
		}

		// Each seed location corresponding to the current K-mer must be
		// used to attempt to extend a match.
		// for seedLoc := seeds; seedLoc != nil; seedLoc = seedLoc.Next { 
		for _, seedLoc := range seeds {
			corSeqId := seedLoc[0]
			corResInd := seedLoc[1]
			corSeq := coarsedb.CoarseSeqGet(corSeqId)

			extCorStart := corResInd + mapSeedSize
			extOrgStart := current + mapSeedSize
			if extCorStart+extSeedSize > corSeq.Len() {
				continue
			}

			cseqExt = corSeq.Residues[extCorStart : extCorStart+extSeedSize]
			oseqExt = orgSeq.Residues[extOrgStart : extOrgStart+extSeedSize]
			if !bytes.Equal(cseqExt, oseqExt) {
				continue
			}

			// The "match" between reference and original sequence will
			// occur somewhere between the the residue index of the seed and
			// the end of the sequence for the reference sequence, and the 
			// position of the "current" pointer and the end of the sequence
			// for the original sequence.
			corMatch, orgMatch := extendMatch(
				corSeq.Residues[corResInd:], orgSeq.Residues[current:],
				db.GappedWindowSize, db.UngappedWindowSize,
				db.MatchKmerSize, db.ExtSeqIdThreshold,
				mem)

			// If the part of the original sequence does not exceed the
			// minimum match length, then we don't accept the match and move
			// on to the next one.
			hasEnd := len(orgMatch)+db.MatchExtend >= orgSeq.Len()-current
			if len(orgMatch) < db.MinMatchLen && !hasEnd {
				continue
			}

			alignment := nwAlign(corMatch, orgMatch, mem)
			id := cablastp.SeqIdentity(alignment[0], alignment[1])
			if id < db.MatchSeqIdThreshold {
				continue
			}

			// If we're close to the end of the original sequence, extend
			// the match to the end.
			changed := false
			if hasEnd {
				orgMatch = orgSeq.Residues[current:]
				changed = true
			}

			// And if we're close to the end of the last match, extend this
			// match backwards.
			if current-lastMatch <= db.MatchExtend {
				orgMatch = orgSeq.Residues[lastMatch : current+len(orgMatch)]
				current = lastMatch
				changed = true
			}

			// If we've extended our match, we need another alignment.
			if changed {
				alignment = nwAlign(corMatch, orgMatch, mem)
			}

			// Otherwise, we accept the first valid match and move on to the 
			// next kmer after the match ends.
			corStart := int(corResInd)
			corEnd := corStart + len(corMatch)
			orgStart := current
			orgEnd := orgStart + len(orgMatch)

			// If there are residues between the end of the last match
			// and the start of this match, then that means no good match
			// could be found for those residues. Thus, they are added to
			// the reference database. (A pathological LinkToReference is
			// created with an empty diff script that points to the added
			// region in the reference database in its entirety.)
			if orgStart-lastMatch > 0 {
				orgSub := orgSeq.NewSubSequence(lastMatch, current)
				orgSubCpy := make([]byte, len(orgSub.Residues))
				copy(orgSubCpy, orgSub.Residues)

				nextCorSeqId := addWithoutMatch(coarsedb, orgSeqId, orgSubCpy)
				cseq.Add(cablastp.NewLinkToCoarseNoDiff(
					nextCorSeqId, 0, len(orgSubCpy)))
			}

			// For the given match, add a LinkToReference to the portion of 
			// the reference sequence matched. This serves as a component 
			// of a compressed original sequence. Also, add a 
			// LinkToCompressed to the reference sequence matched. This 
			// serves as a bridge to expand coarse sequences into their 
			// original sequences.
			cseq.Add(cablastp.NewLinkToCoarse(
				corSeqId, corStart, corEnd, alignment))
			corSeq.AddLink(cablastp.NewLinkToCompressed(
				int32(orgSeqId), int16(corStart), int16(corEnd)))

			// Skip the current pointer ahead to the end of this match.
			// Update the lastMatch pointer to point at the end of this 
			// match.
			lastMatch = orgEnd
			current = orgEnd - 1

			break
		}
	}

	// If there are any leftover residues, then no good match for them
	// could be found. Therefore, add them to the reference database and
	// create the appropriate links.
	if orgSeq.Len()-lastMatch > 0 {
		orgSub := orgSeq.NewSubSequence(lastMatch, orgSeq.Len())
		orgSubCpy := make([]byte, len(orgSub.Residues))
		copy(orgSubCpy, orgSub.Residues)

		nextCorSeqId := addWithoutMatch(coarsedb, orgSeqId, orgSubCpy)
		cseq.Add(cablastp.NewLinkToCoarseNoDiff(
			nextCorSeqId, 0, len(orgSubCpy)))
	}

	return cseq
}

// extendMatch uses a combination of ungapped and gapped extension to find
// quality candidates for compression.
//
// More details to come soon.
func extendMatch(corRes, orgRes []byte,
	gappedWindowSize, ungappedWindowSize, kmerSize, idThreshold int,
	mem *memory) (corMatchRes, orgMatchRes []byte) {

	// Starting at seedLoc.resInd and current, refMatchLen and 
	// orgMatchLen correspond to the length of the match each of
	// the reference and the original sequence, respectively.
	// At the end of the loop, the slices [seedLoc.resInd:refMatchLen]
	// and [current:orgMatchLen] will correspond to the match.
	//
	// A LinkToCoarse is then created with the reference sequence
	// id (seedLoc.seqInd), the start and stop residues of the match
	// in the reference sequence (seedLoc.resInd and
	// (seedLoc.ResInd + refMatchLen)), and an edit script created
	// by an alignment between the matched regions of the reference
	// and original sequences. (A LinkToCoarse corresponds to a
	// single component of a CompressedSeq.)
	corMatchLen, orgMatchLen := 0, 0
	for {
		// If the match has consumed either of the reference or original
		// sequence, then we must quit with what we have.
		if corMatchLen == len(corRes) || orgMatchLen == len(orgRes) {
			break
		}

		// Ungapped extension returns an integer corresponding to the
		// number of residues that the match was extended by.
		matchLen := alignUngapped(
			corRes[corMatchLen:], orgRes[orgMatchLen:],
			ungappedWindowSize, kmerSize, idThreshold)

		// Since ungapped extension increases the reference and
		// original sequence match portions equivalently, add the
		// match length to both.
		corMatchLen += matchLen
		orgMatchLen += matchLen

		// Gapped extension returns an alignment corresponding to the
		// window starting after the previous ungapped extension
		// ended plus the gapped window size. (It is bounded by the
		// length of each sequence.)
		alignment := nwAlign(
			corRes[corMatchLen:min(len(corRes), corMatchLen+gappedWindowSize)],
			orgRes[orgMatchLen:min(len(orgRes), orgMatchLen+gappedWindowSize)],
			mem)

		// If the alignment has a sequence identity below the
		// threshold, then gapped extension has failed. We therefore
		// quit and are forced to be satisfied with whatever
		// refMatchLen and orgMatchLen are set to.
		id := cablastp.SeqIdentity(alignment[0], alignment[1])
		if id < idThreshold {
			break
		}

		// We live to die another day.
		// We need to add to the refMatch{Pos,Len} and orgMatch{Pos,Len}
		// just like we did for ungapped extension. However, an
		// alignment can correspond to two different sized subsequences
		// of the reference and original sequence. Therefore, only
		// increase each by the corresponding sizes from the
		// alignment.
		corMatchLen += alignLen(alignment[0])
		orgMatchLen += alignLen(alignment[1])
	}

	return corRes[:corMatchLen], orgRes[:orgMatchLen]
}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the reference database to the reference database.
// A LinkToCompressed is created and automatically added to the new reference
// sequence.
//
// The id of the new reference sequence added is returned.
func addWithoutMatch(coarsedb *cablastp.CoarseDB, orgSeqId int,
	oseq []byte) int {

	corSeqId, corSeq := coarsedb.Add(oseq)
	corSeq.AddLink(
		cablastp.NewLinkToCompressed(int32(orgSeqId), 0, int16(len(oseq))))
	return corSeqId
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
