package cablastp

import (
	"fmt"
	"log"
	"strings"
	"sync"

	"code.google.com/p/biogo/seq"
)

// SeqIdentity computes the sequence identity of two byte slices.
// The number returned is an integer in the range 0-100, inclusive.
// SeqIdentity returns zero if the lengths of both seq1 and seq2 are zero.
//
// If the lengths of seq1 and seq2 are not equal, SeqIdentity will panic.
func SeqIdentity(seq1, seq2 []byte) int {
	if len(seq1) != len(seq2) {
		log.Panicf("Sequence identity requires that len(seq1) == len(seq2), "+
			"but %d != %d.", len(seq1), len(seq2))
	}
	if len(seq1) == 0 && len(seq2) == 0 {
		return 0
	}

	same := 0
	for i, r1 := range seq1 {
		if r1 == seq2[i] {
			same++
		}
	}
	return (same * 100) / len(seq1)
}

// sequence is the underlying (i.e., embedded) type of reference and original 
// sequences used in cablast.
type sequence struct {
	Name     string
	Residues []byte
	Offset   int
	Id       int
}

// newSeq creates a new sequence and upper cases the given residues.
func newSeq(id int, name string, residues []byte) *sequence {
	residuesStr := strings.ToUpper(string(residues))
	residuesStr = strings.Replace(residuesStr, "*", "", -1)
	return &sequence{
		Name:     name,
		Residues: []byte(residuesStr),
		Offset:   0,
		Id:       id,
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
	s := newSeq(seq.Id, seq.Name, seq.Residues[start:end])
	s.Offset += start
	return s
}

// BiogoSeq returns a new *seq.Seq from biogo.
func (s *sequence) BiogoSeq() *seq.Seq {
	return seq.New(s.Name, s.Residues, nil)
}

// Len retuns the number of residues in this sequence.
func (seq *sequence) Len() int {
	return len(seq.Residues)
}

// String returns a string (fasta) representation of this sequence. If this 
// sequence is a subsequence, then the range of the subsequence (with respect 
// to the original sequence) is also printed.
func (seq *sequence) String() string {
	if seq.Offset == 0 {
		return fmt.Sprintf("> %s (%d)\n%s",
			seq.Name, seq.Id, string(seq.Residues))
	}
	return fmt.Sprintf("> %s (%d) (%d, %d)\n%s",
		seq.Name, seq.Id, seq.Offset, seq.Len(), string(seq.Residues))
}

// referenceSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type ReferenceSeq struct {
	*sequence
	Links    *LinkToCompressed
	linkLock *sync.RWMutex
}

func NewReferenceSeq(id int, name string, residues []byte) *ReferenceSeq {
	return &ReferenceSeq{
		sequence: newSeq(id, name, residues),
		Links:    nil,
		linkLock: &sync.RWMutex{},
	}
}

func NewBiogoReferenceSeq(id int, seq *seq.Seq) *ReferenceSeq {
	return NewReferenceSeq(id, seq.ID, seq.Seq)
}

func (rseq *ReferenceSeq) NewSubSequence(start, end int) *ReferenceSeq {
	return &ReferenceSeq{
		sequence: rseq.sequence.newSubSequence(start, end),
		Links:    nil,
	}
}

func (rseq *ReferenceSeq) AddLink(lk *LinkToCompressed) {
	rseq.linkLock.Lock()
	if rseq.Links == nil {
		rseq.Links = lk
	} else {
		rseq.Links.Next = lk
	}
	rseq.linkLock.Unlock()
}

// OriginalSeq embeds a sequence and serves as a typing mechanism to
// distguish reference sequences in the compressed database with original
// sequences from the input FASTA file.
type OriginalSeq struct {
	*sequence
}

func NewOriginalSeq(id int, name string, residues []byte) *OriginalSeq {
	return &OriginalSeq{sequence: newSeq(id, name, residues)}
}

func NewBiogoOriginalSeq(id int, seq *seq.Seq) *OriginalSeq {
	return &OriginalSeq{sequence: newBiogoSeq(id, seq)}
}

func (oseq *OriginalSeq) NewSubSequence(start, end int) *OriginalSeq {
	return &OriginalSeq{oseq.sequence.newSubSequence(start, end)}
}
