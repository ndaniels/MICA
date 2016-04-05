package mica

import (
	"fmt"
	"log"
	"strings"
	"sync"

	"github.com/TuftsBCB/seq"
)

// SeqIdentity computes the Sequence identity of two byte slices.
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

// Sequence is the underlying (i.e., embedded) type of reference and original
// Sequences used in cablast.
type Sequence struct {
	Name     string
	Residues []byte
	Offset   uint
	Id       int
}

// newSeq creates a new Sequence and upper cases the given residues.
func newSeq(id int, name string, residues []byte) *Sequence {
	residuesStr := strings.ToUpper(string(residues))
	residuesStr = strings.Replace(residuesStr, "*", "", -1)
	return &Sequence{
		Name:     name,
		Residues: []byte(residuesStr),
		Offset:   0,
		Id:       id,
	}
}

// newFastaSeq creates a new *Sequence value from seq's Sequence type, and
// ensures that all residues in the Sequence are upper cased.
func newFastaSeq(id int, s seq.Sequence) *Sequence {
	return newSeq(id, s.Name, s.Bytes())
}

// newSubSequence returns a new *Sequence value that corresponds to a
// subSequence of 'Sequence'. 'start' and 'end' specify an inclusive range in
// 'Sequence'. newSubSequence panics if the range is invalid.
func (seq *Sequence) newSubSequence(start, end uint) *Sequence {
	if start < 0 || start >= end || end > uint(seq.Len()) {
		panic(fmt.Sprintf("Invalid sub Sequence (%d, %d) for Sequence "+
			"with length %d.", start, end, seq.Len()))
	}
	s := newSeq(seq.Id, seq.Name, seq.Residues[start:end])
	s.Offset += start
	return s
}

// FastaSeq returns a new seq.Sequence from TuftsBCB/seq.
func (s *Sequence) FastaSeq() seq.Sequence {
	rs := make([]seq.Residue, len(s.Residues))
	for i := range s.Residues {
		rs[i] = seq.Residue(s.Residues[i])
	}
	return seq.Sequence{s.Name, rs}
}

// Len retuns the number of residues in this Sequence.
func (seq *Sequence) Len() int {
	return len(seq.Residues)
}

// String returns a string (fasta) representation of this Sequence. If this
// Sequence is a subSequence, then the range of the subSequence (with respect
// to the original Sequence) is also printed.
func (seq *Sequence) String() string {
	if seq.Offset == 0 {
		return fmt.Sprintf("> %s (%d)\n%s",
			seq.Name, seq.Id, string(seq.Residues))
	}
	return fmt.Sprintf("> %s (%d) (%d, %d)\n%s",
		seq.Name, seq.Id, seq.Offset, seq.Len(), string(seq.Residues))
}

// referenceSeq embeds a Sequence and serves as a typing mechanism to
// distguish reference Sequences in the compressed database with original
// Sequences from the input FASTA file.
type CoarseSeq struct {
	*Sequence
	Links    *LinkToCompressed
	linkLock *sync.RWMutex
}

func NewCoarseSeq(id int, name string, residues []byte) *CoarseSeq {
	return &CoarseSeq{
		Sequence: newSeq(id, name, residues),
		Links:    nil,
		linkLock: &sync.RWMutex{},
	}
}

func NewFastaCoarseSeq(id int, s seq.Sequence) *CoarseSeq {
	return NewCoarseSeq(id, s.Name, s.Bytes())
}

func (rseq *CoarseSeq) NewSubSequence(start, end uint) *CoarseSeq {
	return &CoarseSeq{
		Sequence: rseq.Sequence.newSubSequence(start, end),
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

// OriginalSeq embeds a Sequence and serves as a typing mechanism to
// distguish reference Sequences in the compressed database with original
// Sequences from the input FASTA file.
type OriginalSeq struct {
	*Sequence
}

func NewOriginalSeq(id int, name string, residues []byte) *OriginalSeq {
	return &OriginalSeq{Sequence: newSeq(id, name, residues)}
}

func NewFastaOriginalSeq(id int, s seq.Sequence) *OriginalSeq {
	return &OriginalSeq{Sequence: newFastaSeq(id, s)}
}

func (oseq *OriginalSeq) NewSubSequence(start, end uint) *OriginalSeq {
	return &OriginalSeq{oseq.Sequence.newSubSequence(start, end)}
}

// ReducedSeq embeds a Sequence and serves as a typing mechanism to
// distguish reduced-alphabet (DNA) Sequences from amino acid Sequences.
type ReducedSeq struct {
	*Sequence
}

//
func NewReducedSeq(oseq *OriginalSeq) *ReducedSeq {
	return &ReducedSeq{Sequence: newSeq(oseq.Sequence.Id,
		oseq.Sequence.Name,
		Reduce(oseq.Sequence.Residues))}
}

func (rseq *ReducedSeq) NewSubSequence(start, end uint) *ReducedSeq {
	return &ReducedSeq{rseq.Sequence.newSubSequence(start, end)}
}