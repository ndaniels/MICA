package main

import (
	"github.com/BurntSushi/cablastp"
)

// An original sequence may result in a combination of the following things:
// 1) Multiple additions to the reference database (multiple reference
// sequences).
// 2) The seeds table updated with any additions to the reference database.
// 3) Multiple LinkToCompressed added to reference sequence link lists.
func compress(refdb *referenceDB, orgSeqId int,
	orgSeq *cablastp.OriginalSequence) *cablastp.CompressedSeq {

	// Keep track of two pointers. 'current' refers to the residue index in the
	// original sequence that extension is currently originating from.
	// 'lastMatch' refers to the residue index of the *end* of the last match
	// with a reference sequence in the compressed database.
	lastMatch, current := 0, 0

	// Iterate through the original sequence a 'kmer' at a time.
	for current = 0; current < orgSeq.Len()-flagSeedSize; current++ {
		kmer := orgSeq.residues[current : current+flagSeedSize]
		if !allUpperAlpha(kmer) {
			continue
		}

		seeds := refdb.seeds.lookup(kmer)
		possibleMatches := make([]match, 0, len(seeds))

		// Each seed location corresponding to the current K-mer must be
		// used to attempt to extend a match. If there is more than one
		// match, the best one is used.
		// (What does "best" mean? Length for now...)
		for _, seedLoc := range seeds {
			refSeqId := seedLoc.seqInd
			refSeq := cdb.seqs[refSeqId]

			// The "match" between reference and original sequence will
			// occur somewhere between the the residue index of the seed and
			// the end of the sequence for the reference sequence, and the 
			// position of the "current" pointer and the end of the sequence
			// for the original sequence.
			refMatchRes, orgMatchRes := extendMatch(
				refSeq.Residues[seedLoc.resInd:],
				orgSeq.Residues[current:])

			// Use an alignment of the matched subsequences to determine
			// whether the match is worthy (i.e., the length is long enough).
			// An alignment is also used to generate an edit script for a
			// LinkToReference.
			alignment := alignGapped(refMatchRes, orgMatchRes)

			// The match is only good if the alignment length is greater
			// than some user-specified length.
			if alignment.Len() >= flagMinMatchLen {
				fmt.Println("seed to (seed + refMatchLen)",
					seedLoc.resInd, seedLoc.resInd+len(refMatchRes))
				fmt.Println("current to (current + orgMatchLen)",
					current, current+len(orgMatchRes))
				fmt.Println("identity",
					cablastp.SeqIdentity(alignment[0].Seq, alignment[1].Seq))
				fmt.Println("> ", rseq.name)
				fmt.Println(string(refMatchRes))
				fmt.Println("> ", origSeq.name)
				fmt.Println(string(orgMatchRes))
				fmt.Println("")
				fmt.Println(alignment)
				fmt.Println("--------------------------------------------")

				possibleMatches = append(possibleMatches,
					match{
						refSeqId: refSeqId,
						refSeq:   refSeq,
						refStart: seedLoc.resInd,
						refEnd:   seedLoc.resInd + len(refMatchRes),
						orgStart: current,
						orgEnd:   current + len(orgMatchRes),
					})
			}
		}
		if len(possibleMatches) > 0 {
			match := bestMatch(possibleMatches)

			// If there are residues between the end of the last match
			// and the start of this match, then that means no good match
			// could be found for those residues. Thus, they are added to
			// the reference database. (A pathological LinkToReference is
			// created with an empty diff script that points to the added
			// region in the reference database in its entirety.)
			if match.orgStart-lastMatch > 0 {
				orgSub := origSeq.newSubSequence(lastMatch, match.orgStart)
				nextRefSeqId := len(refdb.seqs)
				refdb.add(orgSub)
			}

			bestMatch.rseq.addLink(bestMatch.link)
			current = bestMatch.link.original.origEndRes
			lastMatch = current
		}
	}
	sub := origSeq.newSubSequence(lastMatch, origSeq.Len())
	cdb.addToCompressed(sub)
}

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
		alignment := alignGapped(
			refRes[refMatchLen:min(len(refRes), refMatchLen+winSize)],
			orgRes[orgMatchLen:min(len(orgRes), orgMatchLen+winSize)])

		// If the alignment has a sequence identity below the
		// threshold, then gapped extension has failed. We therefore
		// quit and are forced to be satisfied with whatever
		// refMatchLen and orgMatchLen are set to.
		id := cablastp.SeqIdentity(alignment[0].Seq, alignment[1].Seq)
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
		refMatchLen += alignLen(alignment[0].Seq)
		orgMatchLen += alignLen(alignment[1].Seq)
	}

	return refRes[:refMatchLen], orgRes[:orgMatchLen]
}

func bestMatch(matches []match) match {
	bestMatch := matches[0]
	for _, match := range matches[1:] {
		if bestMatch.Less(match) {
			bestMatch = match
		}
	}
	return bestMatch
}
