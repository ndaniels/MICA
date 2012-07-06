package main

import (
	"log"

	"github.com/kortschak/biogo/seq"
	"github.com/kortschak/biogo/util"
	"github.com/kortschak/biogo/align/nw"
)

var lookUpP util.CTL

func init() {
	m := make(map[int]int)
	for i, v := range blosum62Alphabet {
		m[int(v)] = i
	}
	lookUpP = *util.NewCTL(m)
}

func align(rseq *referenceSeq, oseq *originalSeq) seq.Alignment {
	aligner := &nw.Aligner{
		Matrix: blosum62,
		LookUp: lookUpP,
		GapChar: '-',
	}
	alignment, err := aligner.Align(rseq.BiogoSeq(), oseq.BiogoSeq())
	if err != nil {
		log.Panic(err)
	}
	return alignment
}

