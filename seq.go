package main

import (
	"fmt"
	"log"
	"strings"

	"code.google.com/p/biogo/seq"
)

// identity computes the sequence identity of two byte slices.
// The number returned is an integer in the range 0-100, inclusive.
func identity(seq1, seq2 []byte) int {
	if len(seq1) != len(seq2) {
		log.Panicf("Sequence identity requires that len(seq1) == len(seq2), "+
			"but %d != %d.", len(seq1), len(seq2))
	}

	same := 0
	for i, r1 := range seq1 {
		if r1 == seq2[i] {
			same++
		}
	}
	return (same * 100) / len(seq1)
	// return int(100 * (float64(same) / float64(len(seq1)))) 
}

// sequence is the underlying (i.e., embedded) type of reference and original 
// sequences used in cablast.
type sequence struct {
	name     string
	residues []byte
	offset   int
	original *sequence
	id       int
}

// newSeq creates a new sequence and upper cases the given residues.
func newSeq(id int, name string, residues []byte) *sequence {
	return &sequence{
		name:     name,
		residues: []byte(strings.ToUpper(string(residues))),
		offset:   0,
		original: nil,
		id:       id,
	}
}

// newBiogoSeq creates a new *sequence value from biogo's Seq type, and ensures 
// that all residues in the sequence are upper cased.
func newBiogoSeq(id int, s *seq.Seq) *sequence {
	return newSeq(id, s.ID, s.Seq)
}

// newSubSequence returns a new *sequence value that corresponds to a 
// subsequence of 'sequence'. 'start' and 'end' specify an inclusive range in 
// 'sequence'. newSubSequence panics if the range is invalid.
func (seq *sequence) newSubSequence(start, end int) *sequence {
	if start < 0 || start >= end || end > seq.Len() {
		panic(fmt.Sprintf("Invalid sub sequence (%d, %d) for sequence "+
			"with length %d.", start, end, seq.Len()))
	}
	s := newSeq(seq.id, seq.name, seq.residues[start:end])
	s.offset += start
	if seq.original != nil {
		s.original = seq.original
	} else {
		s.original = seq
	}
	return s
}

// BiogoSeq returns a new *seq.Seq from biogo.
func (s *sequence) BiogoSeq() *seq.Seq {
	return seq.New(s.name, s.residues, nil)
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
		return fmt.Sprintf("> %s (%d)\n%s",
			seq.name, seq.id, string(seq.residues))
	}
	return fmt.Sprintf("> %s (%d) (%d, %d)\n%s",
		seq.name, seq.id, seq.offset, seq.Len(), string(seq.residues))
}

// referenceSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type referenceSeq struct {
	*sequence
	links []*linkEntry
}

func newReferenceSeq(id int, name string, residues []byte) *referenceSeq {
	return &referenceSeq{
		sequence: newSeq(id, name, residues),
		links:    make([]*linkEntry, 0),
	}
}

func newBiogoReferenceSeq(id int, seq *seq.Seq) *referenceSeq {
	return newReferenceSeq(id, seq.ID, seq.Seq)
}

func (rseq *referenceSeq) newSubSequence(start, end int) *referenceSeq {
	return &referenceSeq{
		sequence: rseq.sequence.newSubSequence(start, end),
		links:    nil,
	}
}

func (rseq *referenceSeq) addLink(lkEntry *linkEntry) {
	if rseq.original != nil {
		log.Panicf("Cannot add a link to a subsequence of a " +
			"reference sequence.")
	}
	rseq.links = append(rseq.links, lkEntry)
}

// referenceSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type originalSeq struct {
	*sequence
}

func newOriginalSeq(id int, name string, residues []byte) *originalSeq {
	return &originalSeq{sequence: newSeq(id, name, residues)}
}

func newBiogoOriginalSeq(id int, seq *seq.Seq) *originalSeq {
	return &originalSeq{sequence: newBiogoSeq(id, seq)}
}

func (oseq *originalSeq) newSubSequence(start, end int) *originalSeq {
	return &originalSeq{oseq.sequence.newSubSequence(start, end)}
}
