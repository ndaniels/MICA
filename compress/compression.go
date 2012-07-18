package main

import (
	"github.com/BurntSushi/cablastp"
)

// An original sequence may result in a combination of the following things:
// 1) Multiple additions to the reference database (multiple reference
// sequences).
// 2) The seeds table updated with any additions to the reference database.
// 3) Multiple LinkToCompressed added to reference sequence link lists.
func compress(refdb *referenceDB,
	orgSeq *cablastp.OriginalSequence) *cablastp.CompressedSeq

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

		for _, seedLoc := range seeds {
			refSeq := cdb.seqs[seedLoc.seqInd]
			refRes := refSeq.Residues[seedLoc.resInd:]
			orgRes := orgSeq.Residues[current:]

			refMatchLen, orgMatchLen := 0, 0
			for {
				if matchPos == startOseq.Len() {
					break
				}

				subRseq := startRseq.newSubSequence(matchPos, startRseq.Len())
				subOseq := startOseq.newSubSequence(matchPos, startOseq.Len())

				matchLen := alignUngapped(subRseq, subOseq)
				matchPos += matchLen

				tmpRseq := startRseq.newSubSequence(
					matchPos,
					min(startRseq.Len(), matchPos+flagGappedWindowSize))
				tmpOseq := startOseq.newSubSequence(
					matchPos,
					min(startOseq.Len(), matchPos+flagGappedWindowSize))
				alignment := alignGapped(tmpRseq, tmpOseq)
				id := identity(alignment[0].Seq, alignment[1].Seq)
				if id < flagSeqIdThreshold {
					break
				}

				matchPos += tmpOseq.Len()
			}

			if matchPos-current >= flagMinMatchLen {
				subRseq := rseq.newSubSequence(seedLoc.resInd, matchPos)
				subOseq := origSeq.newSubSequence(current, matchPos)
				alignment := alignGapped(subRseq, subOseq)

				fmt.Println("current to matchPos", current, matchPos)
				fmt.Println("identity",
					identity(alignment[0].Seq, alignment[1].Seq))
				fmt.Println("> ", rseq.name)
				fmt.Println(string(subRseq.residues))
				fmt.Println("> ", origSeq.name)
				fmt.Println(string(subOseq.residues))
				fmt.Println("")
				fmt.Println(alignment)
				fmt.Println("--------------------------------------------")

				link := newLinkEntry(seedLoc.resInd, matchPos, subOseq,
					alignment)
				possibleMatches = append(possibleMatches,
					match{
						rseq: rseq,
						link: link,
					})
			}
		}
		if len(possibleMatches) > 0 {
			bestMatch := possibleMatches[0]
			for _, possibleMatch := range possibleMatches[1:] {
				if bestMatch.Less(possibleMatch) {
					bestMatch = possibleMatch
				}
			}

			if bestMatch.link.original.origStartRes-lastMatch > 0 {
				sub := origSeq.newSubSequence(
					lastMatch, bestMatch.link.original.origStartRes)
				fmt.Println(strings.Repeat("#", 45))
				fmt.Println(sub)
				fmt.Println(strings.Repeat("#", 45))
				cdb.addToCompressed(sub)
			}

			bestMatch.rseq.addLink(bestMatch.link)
			current = bestMatch.link.original.origEndRes
			lastMatch = current
		}
	}
	sub := origSeq.newSubSequence(lastMatch, origSeq.Len())
	fmt.Println(strings.Repeat("#", 45))
	fmt.Println(sub)
	fmt.Println(strings.Repeat("#", 45))
	cdb.addToCompressed(sub)
}
