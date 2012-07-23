package main

/*
This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the biogo
package: code.google.com/p/biogo.

It's mostly copied from its original form, but it is optimized specifically
for cablastp to limit allocations and to absolve the need for the biogo/seq.Seq
type.
*/

import (
	"github.com/BurntSushi/cablastp/blosum"

	"code.google.com/p/biogo/align/nw"
	"code.google.com/p/biogo/util"
)

const initPoolSize = 15

var (
	nwLookUpP util.CTL
	aligner   = &nw.Aligner{
		Matrix:  blosum.Matrix62,
		LookUp:  nwLookUpP,
		GapChar: '-',
	}
)

func init() {
	m := make(map[int]int)
	for i, v := range blosum.Alphabet62 {
		m[int(v)] = i
	}
	nwLookUpP = *util.NewCTL(m)
}

func nwAlign(rseq, oseq []byte, table [][]int) [2][]byte {
	gap := len(aligner.Matrix) - 1
	r, c := len(rseq)+1, len(oseq)+1

	if table == nil || len(table) != r {
		table = make([][]int, r)
		for j := range table {
			table[j] = make([]int, c)
		}
	} else {
		for i := range table {
			table[i][i] = 0
		}
	}

	var sdiag, sup, sleft, rVal, oVal int
	valToCode := aligner.LookUp.ValueToCode
	gapChar := aligner.GapChar
	matrix := aligner.Matrix

	for i := 1; i < r; i++ {
		for j := 1; j < c; j++ {
			rVal, oVal = valToCode[rseq[i-1]], valToCode[oseq[j-1]]
			if rVal < 0 || oVal < 0 {
				continue
			} else {
				sdiag = table[i-1][j-1] + matrix[rVal][oVal]
				sup = table[i-1][j] + matrix[rVal][gap]
				sleft = table[i][j-1] + matrix[gap][oVal]
				switch {
				case sdiag > sup && sdiag > sleft:
					table[i][j] = sdiag
				case sup > sleft:
					table[i][j] = sup
				default:
					table[i][j] = sleft
				}
			}
		}
	}

	refAln := make([]byte, 0, len(rseq)+10)
	orgAln := make([]byte, 0, len(oseq)+10)

	i, j := r-1, c-1
	for i > 0 && j > 0 {
		rVal, oVal = valToCode[rseq[i-1]], valToCode[oseq[j-1]]
		if rVal < 0 || oVal < 0 {
			continue
		} else {
			sdiag = table[i-1][j-1] + matrix[rVal][oVal]
			sup = table[i-1][j] + matrix[gap][oVal]
			sleft = table[i][j-1] + matrix[rVal][gap]
			switch {
			case sdiag > sup && sdiag > sleft:
				i--
				j--
				refAln = append(refAln, rseq[i])
				orgAln = append(orgAln, oseq[j])
			case sup > sleft:
				i--
				refAln = append(refAln, rseq[i])
				orgAln = append(orgAln, gapChar)
			default:
				j--
				refAln = append(refAln, gapChar)
				orgAln = append(orgAln, oseq[j])
			}
		}
	}

	for ; i > 0; i-- {
		refAln = append(refAln, rseq[i-1])
		orgAln = append(orgAln, gapChar)
	}
	for ; j > 0; j-- {
		refAln = append(refAln, gapChar)
		orgAln = append(orgAln, oseq[j-1])
	}

	if len(refAln) == len(orgAln) {
		for i, j := 0, len(refAln)-1; i < j; i, j = i+1, j-1 {
			refAln[i], refAln[j] = refAln[j], refAln[i]
			orgAln[i], orgAln[j] = orgAln[j], orgAln[i]
		}
	} else {
		for i, j := 0, len(refAln)-1; i < j; i, j = i+1, j-1 {
			refAln[i], refAln[j] = refAln[j], refAln[i]
		}
		for i, j := 0, len(orgAln)-1; i < j; i, j = i+1, j-1 {
			orgAln[i], orgAln[j] = orgAln[j], orgAln[i]
		}
	}

	return [2][]byte{refAln, orgAln}
}
