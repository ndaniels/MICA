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

func (sequence *sequence) subSequence(start, end int) *sequence {
	if start < 0 || start >= end || end >= sequence.Len() {
		panic(fmt.Sprintf("Invalid sub sequence (%d, %d) for sequence "+
			"with length %d.", start, end, sequence.Len()))
	}
	newId := fmt.Sprintf("%s (%d - %d)", sequence.ID, start, end)
	return newSeq(seq.New(newId, sequence.Seq.Seq[start:end]))
}

type referenceSeq struct {
	*sequence
}

func newReferenceSeq(seq *seq.Seq) *referenceSeq {
	return &referenceSeq{newSeq(seq)}
}

func (refSeq *referenceSeq) subSequence(start, end int) *referenceSeq {
	return newReferenceSeq(refSeq.sequence.subSequence(start, end))
}

type originalSeq struct {
	seqId int
	start, end int
	*sequence
}

func newOriginalSeq(seqId int, seq *seq.Seq) *originalSeq {
	return &originalSeq{
		seqId: seqId,
		start: 0,
		end: seq.Len(),
		newSeq(seq),
	}
}

func (origSeq *originalSeq) subSequence(start, end int) *originalSeq {
	sliced := newOriginalSeq(origSeq.sequence.subSequence(start, end))
	sliced.start, sliced.end = start, end
	return sliced
}
