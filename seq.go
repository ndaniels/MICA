package main

import (
	"fmt"
	"log"
	"strings"

	"github.com/kortschak/biogo/seq"
)

// identity computes the sequence identity of two byte slices.
// The number returned is an integer in the range 0-100, inclusive.
func identity(seq1, seq2 []byte) int {
	if len(seq1) != len(seq2) {
		log.Panicf("Sequence identity requires that len(seq1) == len(seq2), " +
			"but %d != %d.", len(seq1), len(seq2))
	}

	same := 0
	for _, r1 := range seq1 {
		for _, r2 := range seq2 {
			if r1 == r2 {
				same++
			}
		}
	}
	return (same * 100) / len(seq1)
}

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
	if start < 0 || start >= end || end > sequence.Len() {
		panic(fmt.Sprintf("Invalid sub sequence (%d, %d) for sequence "+
			"with length %d.", start, end, sequence.Len()))
	}
	newId := fmt.Sprintf("%s (%d - %d)", sequence.ID, start, end)
	return newSeq(seq.New(newId, sequence.Seq.Seq[start:end], nil))
}

type referenceSeq struct {
	*sequence
}

func newReferenceSeq(seq *seq.Seq) *referenceSeq {
	return &referenceSeq{newSeq(seq)}
}

func (refSeq *referenceSeq) subSequence(start, end int) *referenceSeq {
	return newReferenceSeq(refSeq.sequence.subSequence(start, end).Seq)
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
		sequence: newSeq(seq),
	}
}

func (origSeq *originalSeq) subSequence(start, end int) *originalSeq {
	sliced := newOriginalSeq(origSeq.seqId,
		origSeq.sequence.subSequence(start, end).Seq)
	sliced.start, sliced.end = start, end
	return sliced
}
