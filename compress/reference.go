package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/BurntSushi/cablastp"
)

// referenceDB represents a set of unique sequences that comprise the "coarse"
// database. Sequences in the ReferenceDB are use to re-create the original
// sequences.
type referenceDB struct {
	seqs  []*cablastp.ReferenceSeq
	seeds seeds
}

// newReferenceDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func newReferenceDB() *referenceDB {
	refdb := &referenceDB{
		seqs:  make([]*cablastp.ReferenceSeq, 0, 1000),
		seeds: newSeeds(),
	}
	return refdb
}

// add takes an original sequence, converts it to a reference sequence, and
// adds it as a new reference sequence to the reference database. Seeds are
// also generated for each K-mer in the sequence. The resulting reference
// sequence is returned.
func (refdb *referenceDB) add(
	orgSeq *cablastp.OriginalSeq) *cablastp.ReferenceSeq {

	nextIndex := len(refdb.seqs)
	refSeq := cablastp.NewReferenceSeq(nextIndex, orgSeq.Name, orgSeq.Residues)
	refdb.seeds.add(nextIndex, refSeq)
	refdb.seqs = append(refdb.seqs, refSeq)
	return refSeq
}

// save will save the reference database as a coarse FASTA file and a binary
// encoding of all reference links.
func (refdb *referenceDB) save(fastaFile, linksFile string) error {
	return nil
}

// savePlain will save the reference database as a coarse FASTA file and a 
// plain text encoding of all reference links.
func (refdb *referenceDB) savePlain(fastaFile, linksFile string) error {
	if err := refdb.saveFasta(fastaFile); err != nil {
		return err
	}

	f, err := os.Create(linksFile)
	if err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(f)
	for i, seq := range refdb.seqs {
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

// saveFasta will save the reference database as a coarse FASTA file.
// This should *only* be used for debugging, since a coarse FASTA file is
// usesless without the reference links.
func (refdb *referenceDB) saveFasta(fastaFile string) error {
	f, err := os.Create(fastaFile)
	if err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(f)
	for i, seq := range refdb.seqs {
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
func (refdb *referenceDB) String() string {
	seqStrs := make([]string, len(refdb.seqs))
	for i, seq := range refdb.seqs {
		seqStrs[i] = seq.String()
	}
	return strings.Join(seqStrs, "\n")
}
