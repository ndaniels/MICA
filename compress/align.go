package main

import (
	"log"

	"github.com/BurntSushi/cablastp/blosum"

	"code.google.com/p/biogo/align/nw"
	"code.google.com/p/biogo/seq"
	"code.google.com/p/biogo/util"
)

var lookUpP util.CTL

func init() {
	m := make(map[int]int)
	for i, v := range blosum.Alphabet62 {
		m[int(v)] = i
	}
	lookUpP = *util.NewCTL(m)
}

// alignGapped takes two byte slices and runs Needleman-Wunsch on them to
// form an alignment.
func alignGapped(rseq []byte, oseq []byte) seq.Alignment {
	aligner := &nw.Aligner{
		Matrix:  blosum.Matrix62,
		LookUp:  lookUpP,
		GapChar: '-',
	}
	alignment, err := aligner.Align(seq.Seq{Seq: rseq}, seq.Seq{Seq: oseq})
	if err != nil {
		log.Panic(err)
	}
	return alignment
}

// alignLen computes the length of a sequence in an alignment.
// (i.e., the number of residues that aren't "-".)
func alignLen(seq []byte) (length int) {
	for _, res := range seq {
		if res != '-' {
			length++
		}
	}
	return
}

// alignUngapped takes a reference and an original sub-sequence
// and a starting offset for each sequence. A length
// corresponding to the number of amino acids scanned by greedily consuming
// successive 4-mer matches in 10-mer windows.
//
// The algorithm works by attempting to find *exact* 4-mer matches between the
// sequences in 10-mer windows. If 10 residues are scanned and no 4-mer match
// is found, the the current value of length is returned (which may be 0).
// If a 4-mer match is found, the current value of length is set to the total
// number of amino acid residues scanned, and a search for the next 4-mer match
// for the next 10-mer window is started.
func alignUngapped(rseq []byte, oseq []byte) int {
	length, scanned, successive := 0, 0, 0
	tryNextWindow := true
	for tryNextWindow {
		tryNextWindow = false
		for i := 0; i < flagUngappedWindowSize; i++ {
			// If we've scanned all residues in one of the sub-sequences, then
			// there is nothing left to do for ungapped extension. Therefore,
			// quit and return the number of residues scanned up until the
			// *last* match.
			if scanned >= len(rseq) || scanned >= len(oseq) {
				break
			}

			if rres[scanned] == ores[scanned] {
				successive++
			} else {
				successive = 0
			}

			scanned++
			if successive == flagMatchKmerSize {
				// Get the residues between matches: i.e., after the last
				// match to the start of this match. But only if there is at
				// least one residue in that range.
				if (scanned-flagMatchKmerSize)-length > 0 {
					refBetween := rseq[length : scanned-flagMatchKmerSize]
					orgBetween := oseq[length : scanned-flagMatchKmerSize]
					id := cablastp.SeqIdentity(refBetween, orgBetween)

					// If the identity is less than the threshold, then this
					// K-mer match is no good. But keep trying until the window
					// is closed. (We "keep trying" by decrementing successive
					// matches by 1.)
					if id < flagSeqIdThreshold {
						successive--
						continue
					}
				}

				// If we're here, then we've found a valid match. Update the
				// length to indicate the number of residues scanned and make
				// sure we try the next Ungapped window.
				length = scanned
				successive = 0
				tryNextWindow = true
				break
			}
		}
	}
	return length
}
