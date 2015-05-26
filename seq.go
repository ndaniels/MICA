package cablastp

import (
	"fmt"
	"log"
	"strings"
	"sync"

	"github.com/TuftsBCB/seq"
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

// IsLowComplexity detects whether the residue at the given offset is in
// a region of low complexity, where low complexity is defined as a window
// where every residue is the same (no variation in composition).
func IsLowComplexity(residues []byte, offset, window int) bool {
	repeats := 1
	last := byte(0)
	start := max(0, offset-window)
	end := min(len(residues), offset+window)
	for i := start; i < end; i++ {
		if residues[i] == last {
			repeats++
			if repeats >= window {
				return true
			}
			continue
		}

		// The last residue isn't the same as this residue, so reset.
		last = residues[i]
		repeats = 1
	}
	return false
}

// repetitive returns true if every byte in `bs` is the same.
func repetitive(bs []byte) bool {
	if len(bs) <= 1 {
		return false
	}
	b1 := bs[0]
	for _, b2 := range bs[1:] {
		if b1 != b2 {
			return false
		}
	}
	return true
}

// sequence is the underlying (i.e., embedded) type of reference and original
// sequences used in cablast.
type sequence struct {
	Name     string
	Residues []byte
	Offset   uint
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

// newFastaSeq creates a new *sequence value from seq's Sequence type, and
// ensures that all residues in the sequence are upper cased.
func newFastaSeq(id int, s seq.Sequence) *sequence {
	return newSeq(id, s.Name, s.Bytes())
}

// newSubSequence returns a new *sequence value that corresponds to a
// subsequence of 'sequence'. 'start' and 'end' specify an inclusive range in
// 'sequence'. newSubSequence panics if the range is invalid.
func (seq *sequence) newSubSequence(start, end uint) *sequence {
	if start < 0 || start >= end || end > uint(seq.Len()) {
		panic(fmt.Sprintf("Invalid sub sequence (%d, %d) for sequence "+
			"with length %d.", start, end, seq.Len()))
	}
	s := newSeq(seq.Id, seq.Name, seq.Residues[start:end])
	s.Offset += start
	return s
}

// FastaSeq returns a new seq.Sequence from TuftsBCB/seq.
func (s *sequence) FastaSeq() seq.Sequence {
	rs := make([]seq.Residue, len(s.Residues))
	for i := range s.Residues {
		rs[i] = seq.Residue(s.Residues[i])
	}
	return seq.Sequence{s.Name, rs}
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
type CoarseSeq struct {
	*sequence
	Links    *LinkToCompressed
	linkLock *sync.RWMutex
}

func NewCoarseSeq(id int, name string, residues []byte) *CoarseSeq {
	return &CoarseSeq{
		sequence: newSeq(id, name, residues),
		Links:    nil,
		linkLock: &sync.RWMutex{},
	}
}

func NewFastaCoarseSeq(id int, s seq.Sequence) *CoarseSeq {
	return NewCoarseSeq(id, s.Name, s.Bytes())
}

func (rseq *CoarseSeq) NewSubSequence(start, end uint) *CoarseSeq {
	return &CoarseSeq{
		sequence: rseq.sequence.newSubSequence(start, end),
		Links:    nil,
	}
}

func (rseq *CoarseSeq) AddLink(link *LinkToCompressed) {
	rseq.linkLock.Lock()
	rseq.addLink(link)
	rseq.linkLock.Unlock()
}

func (rseq *CoarseSeq) addLink(link *LinkToCompressed) {
	if rseq.Links == nil {
		rseq.Links = link
	} else {
		lk := rseq.Links
		for ; lk.Next != nil; lk = lk.Next {
		}
		lk.Next = link
	}
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

func NewFastaOriginalSeq(id int, s seq.Sequence) *OriginalSeq {
	return &OriginalSeq{sequence: newFastaSeq(id, s)}
}

func (oseq *OriginalSeq) NewSubSequence(start, end uint) *OriginalSeq {
	return &OriginalSeq{oseq.sequence.newSubSequence(start, end)}
}
