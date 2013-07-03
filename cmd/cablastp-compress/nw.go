package main

import (
	"github.com/BurntSushi/cablastp/blosum"
)

var (
	resTrans [256]int
)

// Initialize the alignment lookup table. (i.e., translate ASCII residue
// characters to BLOSUM62 matrix indices.)
func init() {
	for i := 0; i < len(blosum.Alphabet62); i++ {
		resTrans[blosum.Alphabet62[i]] = i
	}
}

// nwAlign performs Needleman-Wunsch sequence alignment.
//
// This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the
// biogo package: code.google.com/p/biogo.
//
// It's mostly copied from its original form, but it is optimized specifically
// for cablastp to limit allocations and to absolve the need for the
// biogo/seq.Seq type. There are several additional optimizations to limit
// functional calls and pointer deferences.
//
// Perhaps the biggest optimization; however, is constraining dynamic
// programming to only allow a limited number of gaps proportion to the
// length of the large of rseq and oseq.
func nwAlign(rseq, oseq []byte, mem *memory) [2][]byte {
	gap := len(blosum.Matrix62) - 1
	r, c := len(rseq)+1, len(oseq)+1
	off := 0

	constrained := true
	constraint := max(r, c) / 4
	if r <= 11 || c <= 11 {
		constrained = false
	}

	var table []int
	var j int
	if r*c > dynamicTableSize {
		table = make([]int, r*c)
	} else {
		table = mem.table[:r*c]
		for i := range table {
			table[i] = 0
		}
	}

	var sdiag, sup, sleft, rVal, oVal int
	gapChar := byte('-')
	matrix := blosum.Matrix62

	var i2, i3 int
	for i := 1; i < r; i++ {
		i2 = (i - 1) * c
		i3 = i * c
		for j = 1; j < c; j++ {
			if constrained && ((i-j) > constraint || (j-i) > constraint) {
				continue
			}
			rVal, oVal = resTrans[rseq[i-1]], resTrans[oseq[j-1]]

			off = i2 + (j - 1)
			sdiag = table[off] + matrix[rVal][oVal]
			sup = table[off+1] + matrix[rVal][gap]
			sleft = table[off+c] + matrix[gap][oVal]
			switch {
			case sdiag > sup && sdiag > sleft:
				table[i3+j] = sdiag
			case sup > sleft:
				table[i3+j] = sup
			default:
				table[i3+j] = sleft
			}
		}
	}

	refAln, orgAln := mem.ref[:0], mem.org[:0]

	i, j := r-1, c-1
	for i > 0 && j > 0 {
		rVal, oVal = resTrans[rseq[i-1]], resTrans[oseq[j-1]]

		sdiag = table[(i-1)*c+(j-1)] + matrix[rVal][oVal]
		sup = table[(i-1)*c+j] + matrix[gap][oVal]
		sleft = table[i*c+(j-1)] + matrix[rVal][gap]
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
