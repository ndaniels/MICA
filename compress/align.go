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

func alignGapped(rseq *referenceSeq, oseq *originalSeq) seq.Alignment {
	aligner := &nw.Aligner{
		Matrix:  blosum.Matrix62,
		LookUp:  lookUpP,
		GapChar: '-',
	}
	alignment, err := aligner.Align(rseq.BiogoSeq(), oseq.BiogoSeq())
	if err != nil {
		log.Panic(err)
	}
	return alignment
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
func alignUngapped(rseq *referenceSeq, oseq *originalSeq) int {
	rres, ores := rseq.residues, oseq.residues
	length, scanned, successive := 0, 0, 0
	tryNext10 := true
	for tryNext10 {
		tryNext10 = false
		for i := 0; i < 10; i++ {
			if scanned >= len(rres) || scanned >= len(ores) {
				break
			}
			if rres[scanned] == ores[scanned] {
				successive++
			} else {
				successive = 0
			}

			scanned++
			if successive == flagMatchKmerSize {
				beforeMatch := scanned - flagMatchKmerSize
				id := identity(rres[length:beforeMatch],
					ores[length:beforeMatch])
				if id < flagSeqIdThreshold {
					break
				}
				length = scanned
				successive = 0
				tryNext10 = true
				break
			}
		}
	}
	return length
}
