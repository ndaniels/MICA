package main

import (
	"github.com/kortschak/biogo/util"
	"github.com/kortschak/biogo/align/sw"
)

var lookUpP util.CTL

func init() {
	m := make(map[int]int)
	for i, v := range blosum62Alphabet {
		m[int(v)] = i
	}
	lookUpP = *util.NewCTL(m)
}

func align(rseq *referenceSeq, oseq *referenceSeq) seq.Alignment {
}

