package cablastp

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// CoarseDB represents a set of unique sequences that comprise the "coarse"
// database. Sequences in the ReferenceDB are use to re-create the original
// sequences.
type CoarseDB struct {
	Seqs  []*ReferenceSeq
	Seeds *Seeds
}

// NewCoarseDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func NewCoarseDB(seedSize int) *CoarseDB {
	coarsedb := &CoarseDB{
		Seqs:  make([]*ReferenceSeq, 0, 1000),
		Seeds: NewSeeds(seedSize),
	}
	return coarsedb
}

// Add takes an original sequence, converts it to a reference sequence, and
// adds it as a new reference sequence to the reference database. Seeds are
// also generated for each K-mer in the sequence. The resulting reference
// sequence is returned.
func (coarsedb *CoarseDB) Add(orgSeq *OriginalSeq) *ReferenceSeq {

	nextIndex := len(coarsedb.Seqs)
	refSeq := NewReferenceSeq(nextIndex, orgSeq.Name, orgSeq.Residues)
	coarsedb.Seeds.Add(nextIndex, refSeq)
	coarsedb.Seqs = append(coarsedb.Seqs, refSeq)
	return refSeq
}

// Save will save the reference database as a coarse FASTA file and a binary
// encoding of all reference links.
func (coarsedb *CoarseDB) Save(fastaFile, linksFile string) error {
	return nil
}

// SavePlain will save the reference database as a coarse FASTA file and a 
// plain text encoding of all reference links.
func (coarsedb *CoarseDB) SavePlain(fastaFile, linksFile string) error {
	if err := coarsedb.SaveFasta(fastaFile); err != nil {
		return err
	}

	f, err := os.Create(linksFile)
	if err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(f)
	for i, seq := range coarsedb.Seqs {
		_, err = fmt.Fprintf(bufWriter, "> %d\n", i)
		if err != nil {
			return err
		}
		for _, link := range seq.Links {
			_, err := fmt.Fprintf(bufWriter, "%s\n", link)
			if err != nil {
				return err
			}
		}
	}
	if err := bufWriter.Flush(); err != nil {
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}
	return nil
}

// SaveFasta will save the reference database as a coarse FASTA file.
// This should *only* be used for debugging, since a coarse FASTA file is
// usesless without the reference links.
func (coarsedb *CoarseDB) SaveFasta(fastaFile string) error {
	f, err := os.Create(fastaFile)
	if err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(f)
	for i, seq := range coarsedb.Seqs {
		_, err = fmt.Fprintf(bufWriter,
			"> %d\n%s\n", i, string(seq.Residues))
		if err != nil {
			return err
		}
	}
	if err := bufWriter.Flush(); err != nil {
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}
	return nil
}

// String returns a FASTA representation of the reference database.
func (coarsedb *CoarseDB) String() string {
	seqStrs := make([]string, len(coarsedb.Seqs))
	for i, seq := range coarsedb.Seqs {
		seqStrs[i] = seq.String()
	}
	return strings.Join(seqStrs, "\n")
}
