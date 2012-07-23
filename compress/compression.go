package main

import (
	"sync"

	"github.com/BurntSushi/cablastp"
)

type compressJob struct {
	orgSeqId int
	orgSeq   *cablastp.OriginalSeq
}

func compressWorker(DB *cablastp.DB, jobs chan compressJob,
	wg *sync.WaitGroup) {

	for job := range jobs {
		comSeq := compress(DB.CoarseDB, job.orgSeqId, job.orgSeq)
		DB.ComDB.Write(comSeq)
	}
	wg.Done()
}

func compress(coarsedb *cablastp.CoarseDB, orgSeqId int,
	orgSeq *cablastp.OriginalSeq) cablastp.CompressedSeq {

	cseq := cablastp.NewCompressedSeq(orgSeqId, orgSeq.Name)
	seedSize := coarsedb.Seeds.SeedSize

	// Keep track of two pointers. 'current' refers to the residue index in the
	// original sequence that extension is currently originating from.
	// 'lastMatch' refers to the residue index of the *end* of the last match
	// with a reference sequence in the compressed database.
	lastMatch, current := 0, 0

	// Iterate through the original sequence a 'kmer' at a time.
	for current = 0; current < orgSeq.Len()-seedSize; current++ {
		kmer := orgSeq.Residues[current : current+seedSize]
		if !cablastp.KmerAllUpperAlpha(kmer) {
			continue
		}

		seeds := coarsedb.Seeds.Lookup(kmer)

		// Each seed location corresponding to the current K-mer must be
		// used to attempt to extend a match.
		for _, seedLoc := range seeds {
			refSeqId := seedLoc[0]
			refSeq := coarsedb.RefSeqGet(refSeqId)

			// The "match" between reference and original sequence will
			// occur somewhere between the the residue index of the seed and
			// the end of the sequence for the reference sequence, and the 
			// position of the "current" pointer and the end of the sequence
			// for the original sequence.
			refMatch, orgMatch := extendMatch(
				refSeq.Residues[seedLoc[1]:],
				orgSeq.Residues[current:])

			// If the part of the original sequence does not exceed the
			// minimum match length, then we don't accept the match and move
			// on to the next one.
			if len(orgMatch) < flagMinMatchLen {
				continue
			}

			// Otherwise, we accept the first valid match and move on to the 
			// next kmer after the match ends.
			refStart := seedLoc[1]
			refEnd := refStart + len(refMatch)
			orgStart := current
			orgEnd := orgStart + len(orgMatch)
			alignment := nwAlign(refMatch, orgMatch)

			// If there are residues between the end of the last match
			// and the start of this match, then that means no good match
			// could be found for those residues. Thus, they are added to
			// the reference database. (A pathological LinkToReference is
			// created with an empty diff script that points to the added
			// region in the reference database in its entirety.)
			if orgStart-lastMatch > 0 {
				orgSub := orgSeq.NewSubSequence(lastMatch, current)
				nextRefSeqId := addWithoutMatch(coarsedb, orgSeqId, orgSub)
				cseq.Add(cablastp.NewLinkToReferenceNoDiff(
					nextRefSeqId, 0, orgSub.Len()))
			}

			// For the given match, add a LinkToReference to the portion of 
			// the reference sequence matched. This serves as a component 
			// of a compressed original sequence. Also, add a 
			// LinkToCompressed to the reference sequence matched. This 
			// serves as a bridge to expand coarse sequences into their 
			// original sequences.
			cseq.Add(cablastp.NewLinkToReference(
				refSeqId, refStart, refEnd, alignment))
			refSeq.AddLink(cablastp.NewLinkToCompressed(
				orgSeqId, refStart, refEnd))

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
		nextRefSeqId := addWithoutMatch(coarsedb, orgSeqId, orgSub)
		cseq.Add(cablastp.NewLinkToReferenceNoDiff(
			nextRefSeqId, 0, orgSub.Len()))
	}

	return cseq
}

// extendMatch uses a combination of ungapped and gapped extension to find
// quality candidates for compression.
//
// More details to come soon.
func extendMatch(refRes, orgRes []byte) (refMatchRes, orgMatchRes []byte) {
	// Starting at seedLoc.resInd and current, refMatchLen and 
	// orgMatchLen correspond to the length of the match each of
	// the reference and the original sequence, respectively.
	// At the end of the loop, the slices [seedLoc.resInd:refMatchLen]
	// and [current:orgMatchLen] will correspond to the match.
	//
	// A LinkToReference is then created with the reference sequence
	// id (seedLoc.seqInd), the start and stop residues of the match
	// in the reference sequence (seedLoc.resInd and
	// (seedLoc.ResInd + refMatchLen)), and an edit script created
	// by an alignment between the matched regions of the reference
	// and original sequences. (A LinkToReference corresponds to a
	// single component of a CompressedSeq.)
	refMatchLen, orgMatchLen := 0, 0
	for {
		// If the match has consumed either of the reference or original
		// sequence, then we must quit with what we have.
		if refMatchLen == len(refRes) || orgMatchLen == len(orgRes) {
			break
		}

		// Ungapped extension returns an integer corresponding to the
		// number of residues that the match was extended by.
		matchLen := alignUngapped(
			refRes[refMatchLen:], orgRes[orgMatchLen:])

		// Since ungapped extension increases the reference and
		// original sequence match portions equivalently, add the
		// match length to both.
		refMatchLen += matchLen
		orgMatchLen += matchLen

		// Gapped extension returns an alignment corresponding to the
		// window starting after the previous ungapped extension
		// ended plus the gapped window size. (It is bounded by the
		// length of each sequence.)
		winSize := flagGappedWindowSize
		alignment := nwAlign(
			refRes[refMatchLen:min(len(refRes), refMatchLen+winSize)],
			orgRes[orgMatchLen:min(len(orgRes), orgMatchLen+winSize)])
		// If the alignment has a sequence identity below the
		// threshold, then gapped extension has failed. We therefore
		// quit and are forced to be satisfied with whatever
		// refMatchLen and orgMatchLen are set to.
		id := cablastp.SeqIdentity(alignment[0], alignment[1])
		if id < flagSeqIdThreshold {
			break
		}

		// We live to die another day.
		// We need to add to the refMatch{Pos,Len} and orgMatch{Pos,Len}
		// just like we did for ungapped extension. However, an
		// alignment can correspond to two different sized subsequences
		// of the reference and original sequence. Therefore, only
		// increase each by the corresponding sizes from the
		// alignment.
		refMatchLen += alignLen(alignment[0])
		orgMatchLen += alignLen(alignment[1])
	}

	return refRes[:refMatchLen], orgRes[:orgMatchLen]
}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the reference database to the reference database.
// A LinkToCompressed is created and automatically added to the new reference
// sequence.
//
// The id of the new reference sequence added is returned.
func addWithoutMatch(coarsedb *cablastp.CoarseDB, orgSeqId int,
	orgSeq *cablastp.OriginalSeq) int {

	refSeqId, refSeq := coarsedb.Add(orgSeq)
	refSeq.AddLink(cablastp.NewLinkToCompressed(orgSeqId, 0, orgSeq.Len()))
	return refSeqId
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
