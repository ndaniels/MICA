package main

import (
	"strings"

	"github.com/kortschak/biogo/seq"
)

type sequence struct {
	*seq.Seq
}

func newSeq(s *seq.Seq) *sequence {
	s.Seq = []byte(strings.ToUpper(string(s.Seq)))
	return &sequence{s}
}

func (sequence *sequence) residues() []byte {
	return sequence.Seq.Seq
}

type referenceSeq struct {
	*sequence
}

func newReferenceSeq(seq *seq.Seq) *referenceSeq {
	return &referenceSeq{newSeq(seq)}
}

type originalSeq struct {
	*sequence
}

func newOriginalSeq(seq *seq.Seq) *originalSeq {
	return &originalSeq{newSeq(seq)}
}
