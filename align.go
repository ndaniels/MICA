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

func align(rseq *referenceSeq, oseq *originalSeq) seq.Alignment {
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
