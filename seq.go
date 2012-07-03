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

// sequence is the underlying (i.e., embedded) type of reference and original 
// sequences used in cablast.
type sequence struct {
	name string
	residues []byte
	offset int
	original *sequence
}

// newSeq creates a new sequence and upper cases the given residues.
func newSeq(name string, residues []byte) *sequence {
	return &sequence{
		name: name,
		residues: []byte(strings.ToUpper(string(residues))),
		offset: 0,
		original: nil,
	}
}

// newBiogoSeq creates a new *sequence value from biogo's Seq type, and ensures 
// that all residues in the sequence are upper cased.
func newBiogoSeq(s *seq.Seq) *sequence {
	return newSeq(s.ID, s.Seq)
}

// newSubSequence returns a new *sequence value that corresponds to a 
// subsequence of 'sequence'. 'start' and 'end' specify an inclusive range in 
// 'sequence'. newSubSequence panics if the range is invalid.
func (seq *sequence) newSubSequence(start, end int) *sequence {
	if start < 0 || start >= end || end > seq.Len() {
		panic(fmt.Sprintf("Invalid sub sequence (%d, %d) for sequence "+
			"with length %d.", start, end, seq.Len()))
	}
	s := newSeq(seq.Name, seq.residues[start:end])
	s.offset += start
	if seq.original != nil {
		s.original = seq.original
	} else {
		s.original = seq
	}
	return s
}

// Len retuns the number of residues in this sequence.
func (seq *sequence) Len() int {
	return len(seq.residues)
}

// String returns a string (fasta) representation of this sequence. If this 
// sequence is a subsequence, then the range of the subsequence (with respect 
// to the original sequence) is also printed.
func (seq *sequence) String() string {
	if seq.offset == 0 {
		return fmt.Sprintf("> %s\n%s", seq.name, string(seq.residues))
	}
	return fmt.Sprintf("> %s (%d, %d)\n%s",
		seq.name, seq.offset, seq.Len(), string(seq.residues))
}

// referenceSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type referenceSeq struct {
	*sequence
}

func newReferenceSeq(name string, residues []byte) *referenceSeq {
	return &referenceSeq{newSeq(name, residues)}

func newBiogoReferenceSeq(seq *seq.Seq) *referenceSeq {
	return &referenceSeq{newBiogoSeq(seq)}
}

func (rseq *referenceSeq) newSubSequence(start, end int) *referenceSeq {
	return &referenceSeq{rseq.sequence.newSubSequence(start, end)}
}

// referenceSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type originalSeq struct {
	*sequence
}

func newOriginalSeq(name string, residues []byte) *originalSeq {
	return &originalSeq{newSeq(name, residues)}

func newBiogoOriginalSeq(seq *seq.Seq) *originalSeq {
	return &originalSeq{newBiogoSeq(seq)}
}

func (oseq *originalSeq) newSubSequence(start, end int) *originalSeq {
	return &originalSeq{oseq.sequence.newSubSequence(start, end)}
}
