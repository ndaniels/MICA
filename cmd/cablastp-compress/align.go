package main

import (
	"github.com/BurntSushi/cablastp"
)

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

// alignUngapped takes a coarse and an original sub-sequence and returns a
// length corresponding to the number of amino acids scanned by greedily
// consuming successive K-mer matches in N-mer windows.
//
// The algorithm works by attempting to find *exact* K-mer matches between the
// sequences in N-mer windows. If N residues are scanned and no K-mer match
// is found, the the current value of length is returned (which may be 0).
// If a K-mer match is found, the current value of length is set to the total
// number of amino acid residues scanned, and a search for the next K-mer match
// for the next N-mer window is started.
func alignUngapped(rseq []byte, oseq []byte,
	windowSize, kmerSize, idThreshold int) int {

	length, scanned, successive := 0, 0, 0
	tryNextWindow := true
	for tryNextWindow {
		tryNextWindow = false
		for i := 0; i < windowSize; i++ {
			// If we've scanned all residues in one of the sub-sequences, then
			// there is nothing left to do for ungapped extension. Therefore,
			// quit and return the number of residues scanned up until the
			// *last* match.
			if scanned >= len(rseq) || scanned >= len(oseq) {
				break
			}

			if rseq[scanned] == oseq[scanned] {
				successive++
			} else {
				successive = 0
			}

			scanned++
			if successive == kmerSize {
				// Get the residues between matches: i.e., after the last
				// match to the start of this match. But only if there is at
				// least one residue in that range.
				if (scanned-kmerSize)-length > 0 {
					id := cablastp.SeqIdentity(
						rseq[length:scanned-kmerSize],
						oseq[length:scanned-kmerSize])

					// If the identity is less than the threshold, then this
					// K-mer match is no good. But keep trying until the window
					// is closed. (We "keep trying" by decrementing successive
					// matches by 1.)
					if id < idThreshold {
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
