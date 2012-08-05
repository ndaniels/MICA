package cablastp

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"sync"
)

// CoarseDB represents a set of unique sequences that comprise the "coarse"
// database. Sequences in the ReferenceDB are use to re-create the original
// sequences.
type CoarseDB struct {
	Seqs  []*CoarseSeq
	Seeds Seeds

	FileFasta *os.File
	FileSeeds *os.File
	FileLinks *os.File

	seqLock *sync.RWMutex
}

// NewCoarseDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func NewCoarseDB(fastaFile, seedsFile, linksFile *os.File,
	seedSize int) *CoarseDB {

	coarsedb := &CoarseDB{
		Seqs:      make([]*CoarseSeq, 0, 10000000),
		Seeds:     NewSeeds(seedSize),
		FileFasta: fastaFile,
		FileSeeds: seedsFile,
		FileLinks: linksFile,
		seqLock:   &sync.RWMutex{},
	}
	return coarsedb
}

func NewAppendCoarseDB(fastaFile, seedsFile, linksFile *os.File,
	seedSize int) *CoarseDB {

	panic("not implemented")
}

func LoadCoarseDB(fastaFile, linksFile *os.File,
	seedSize int) (*CoarseDB, error) {

	return nil, fmt.Errorf("Not implemented.")
}

// Add takes an original sequence, converts it to a reference sequence, and
// adds it as a new reference sequence to the reference database. Seeds are
// also generated for each K-mer in the sequence. The resulting reference
// sequence is returned.
func (coarsedb *CoarseDB) Add(oseq []byte) (int, *CoarseSeq) {
	coarsedb.seqLock.Lock()
	id := len(coarsedb.Seqs)
	corSeq := NewCoarseSeq(id, "", oseq)
	coarsedb.Seqs = append(coarsedb.Seqs, corSeq)
	coarsedb.seqLock.Unlock()

	coarsedb.Seeds.Add(id, corSeq)

	return id, corSeq
}

func (coarsedb *CoarseDB) CoarseSeqGet(i int) *CoarseSeq {
	coarsedb.seqLock.RLock()
	seq := coarsedb.Seqs[i]
	coarsedb.seqLock.RUnlock()

	return seq
}

func (coarsedb *CoarseDB) ReadClose() {
	coarsedb.FileFasta.Close()
	coarsedb.FileLinks.Close()
}

func (coarsedb *CoarseDB) WriteClose() {
	coarsedb.FileFasta.Close()
	coarsedb.FileSeeds.Close()
	coarsedb.FileLinks.Close()
}

// Save will save the reference database as a coarse FASTA file and a binary
// encoding of all reference links.
func (coarsedb *CoarseDB) Save() error {
	coarsedb.seqLock.RLock()
	defer coarsedb.seqLock.RUnlock()

	return coarsedb.save()
}

// savePlain will save the reference database as a coarse FASTA file and a 
// plain text encoding of all reference links.
func (coarsedb *CoarseDB) save() error {
	if err := coarsedb.saveFasta(); err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(coarsedb.FileLinks)
	for i, seq := range coarsedb.Seqs {
		_, err := fmt.Fprintf(bufWriter, "> %d\n", i)
		if err != nil {
			return err
		}
		for link := seq.Links; link != nil; link = link.Next {
			_, err := fmt.Fprintf(bufWriter, "%s\n", link)
			if err != nil {
				return err
			}
		}
	}
	if err := bufWriter.Flush(); err != nil {
		return err
	}
	return nil
}

// SaveFasta will save the reference database as a coarse FASTA file.
// This should *only* be used for debugging, since a coarse FASTA file is
// usesless without the reference links.
func (coarsedb *CoarseDB) saveFasta() error {
	bufWriter := bufio.NewWriter(coarsedb.FileFasta)
	for i, seq := range coarsedb.Seqs {
		_, err := fmt.Fprintf(bufWriter,
			"> %d\n%s\n", i, string(seq.Residues))
		if err != nil {
			return err
		}
	}
	if err := bufWriter.Flush(); err != nil {
		return err
	}
	return nil
}

// String returns a FASTA representation of the reference database.
func (coarsedb *CoarseDB) String() string {
	coarsedb.seqLock.RLock()
	defer coarsedb.seqLock.RUnlock()

	seqStrs := make([]string, len(coarsedb.Seqs))
	for i, seq := range coarsedb.Seqs {
		seqStrs[i] = seq.String()
	}
	return strings.Join(seqStrs, "\n")
}
