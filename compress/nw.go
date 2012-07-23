package main

/*
This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the biogo
package: code.google.com/p/biogo.

It's mostly copied from its original form, but it is optimized specifically
for cablastp to limit allocations and to absolve the need for the biogo/seq.Seq
type.
*/

import (
	"sync"

	"github.com/BurntSushi/cablastp/blosum"

	"code.google.com/p/biogo/align/nw"
	"code.google.com/p/biogo/util"
)

const (
	nwDiag = iota
	nwUp
	nwLeft
)

const initPoolSize = 15

var (
	nwLookUpP util.CTL
	aligner   = &nw.Aligner{
		Matrix:  blosum.Matrix62,
		LookUp:  nwLookUpP,
		GapChar: '-',
	}

	tablePool = make([][][]int, initPoolSize, 100)
	tableLock = &sync.Mutex{}
)

func init() {
	tableLock.Lock()
	r, c := flagGappedWindowSize, flagGappedWindowSize
	for i := 0; i < initPoolSize; i++ {
		tablePool[i] = make([][]int, r)
		for j := range tablePool[i] {
			tablePool[i][j] = make([]int, c)
		}
	}
	tableLock.Unlock()

	m := make(map[int]int)
	for i, v := range blosum.Alphabet62 {
		m[int(v)] = i
	}
	nwLookUpP = *util.NewCTL(m)
}

func allocTable(r, c int) [][]int {
	if r != flagGappedWindowSize || c != flagGappedWindowSize {
		table := make([][]int, r)
		for j := range table {
			table[j] = make([]int, c)
		}
		return table
	}

	tableLock.Lock()
	for i := len(tablePool) - 1; i >= 0; i-- {
		table := tablePool[i]
		if len(table) == r && len(table[0]) == c {
			tablePool = append(tablePool[:i], tablePool[i+1:]...)
			tableLock.Unlock()
			return table
		}
	}
	tableLock.Unlock()

	// We couldn't find a table with the desired size, so allocate one.
	// It will be added back to the pool when 'freeTable' is called.
	println("Allocating new table:", r, c)
	table := make([][]int, r)
	for j := range table {
		table[j] = make([]int, c)
	}
	return table
}

func freeTable(table [][]int) {
	if len(table) != flagGappedWindowSize ||
		len(table[0]) != flagGappedWindowSize {

		return
	}

	tableLock.Lock()
	defer tableLock.Unlock()

	tablePool = append(tablePool, table)
}

func nwAlign(rseq, oseq []byte) [2][]byte {
	gap := len(aligner.Matrix) - 1
	r, c := len(rseq)+1, len(oseq)+1
	table := allocTable(r, c)

	var sdiag, sup, sleft int
	valToCode := aligner.LookUp.ValueToCode
	gapChar := aligner.GapChar
	matrix := aligner.Matrix

	for i := 1; i < r; i++ {
		for j := 1; j < c; j++ {
			rVal, oVal := valToCode[rseq[i-1]], valToCode[oseq[j-1]]
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

	refAln := make([]byte, 0, len(rseq))
	orgAln := make([]byte, 0, len(oseq))

	i, j := r-1, c-1
	for i > 0 && j > 0 {
		rVal, oVal := valToCode[rseq[i-1]], valToCode[oseq[j-1]]
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

	freeTable(table)

	for ; i > 0; i-- {
		refAln = append(refAln, rseq[i-1])
		orgAln = append(orgAln, gapChar)
	}
	for ; j > 0; j-- {
		refAln = append(refAln, gapChar)
		orgAln = append(orgAln, oseq[j-1])
	}

	for i, j := 0, len(refAln)-1; i < j; i, j = i+1, j-1 {
		refAln[i], refAln[j] = refAln[j], refAln[i]
	}
	for i, j := 0, len(orgAln)-1; i < j; i, j = i+1, j-1 {
		orgAln[i], orgAln[j] = orgAln[j], orgAln[i]
	}

	return [2][]byte{refAln, orgAln}
}
