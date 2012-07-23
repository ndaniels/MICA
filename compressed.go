package cablastp

import (
	"bufio"
	"fmt"
	"os"
)

type CompressedDB struct {
	Seqs       []*CompressedSeq
	File       *os.File
	Index      *os.File
	writerChan chan *CompressedSeq
	writerDone chan struct{}
}

func NewCompressedDB(file *os.File, index *os.File, plain bool) *CompressedDB {
	cdb := &CompressedDB{
		Seqs:       make([]*CompressedSeq, 0, 100),
		File:       file,
		Index:      index,
		writerChan: make(chan *CompressedSeq, 200),
		writerDone: make(chan struct{}, 0),
	}

	if plain {
		go cdb.writerPlain()
	} else {
		go cdb.writerBinary()
	}

	return cdb
}

func (comdb *CompressedDB) Close() {
	close(comdb.writerChan) // will close comdb.File

	// Wait for the writer goroutine to finish.
	<-comdb.writerDone

	comdb.Index.Close()
}

func (comdb *CompressedDB) Add(comSeq *CompressedSeq) {
	comdb.Seqs = append(comdb.Seqs, comSeq)
}

func (comdb *CompressedDB) Len() int {
	return len(comdb.Seqs)
}

func (comdb *CompressedDB) writerPlain() {
	for cseq := range comdb.writerChan {
		_, err := fmt.Fprintf(comdb.File, "> %d; %s\n", cseq.Id, cseq.Name)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}
		for _, link := range cseq.Links {
			_, err := fmt.Fprintf(comdb.File, "%s\n", link)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err)
				os.Exit(1)
			}
		}
	}
	comdb.File.Close()
	comdb.writerDone <- struct{}{}
}

func (comdb *CompressedDB) writerBinary() {
	comdb.File.Close()
}

func (comdb *CompressedDB) Write(cseq *CompressedSeq) {
	comdb.writerChan <- cseq
}

// Save will save the compressed database as a binary encoding of all
// compressed sequences. (A compressed sequence is simply an ordered list of
// links into the reference database.)
func (comdb *CompressedDB) Save(name string) error {
	return nil
}

// SavePlain will save the compressed database as a plain text encoding of all
// compressed sequences. (A compressed sequence is simply an ordered list of
// links into the reference database.)
func (comdb *CompressedDB) SavePlain(name string) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}

	bufWriter := bufio.NewWriter(f)
	for i, seq := range comdb.Seqs {
		_, err = fmt.Fprintf(bufWriter, "> %d; %s\n", i, seq.Name)
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

// CompressedSeq corresponds to the components of a compressed sequence.
type CompressedSeq struct {
	Id int

	// Name is an uncompressed string from the original FASTA header.
	Name string

	// Links is an ordered lists of links to portions of the reference
	// database. When all links are followed, the concatenation of each
	// sequence correspond to each link equals the entire original sequence.
	Links []*LinkToReference
}

// NewCompressedSeq creates a CompressedSeq value using the name provided.
// The Link slice is initialized but empty.
func NewCompressedSeq(id int, name string) *CompressedSeq {
	return &CompressedSeq{
		Id:    id,
		Name:  name,
		Links: make([]*LinkToReference, 0, 15),
	}
}

// Add will add a LinkToReference to the end of the CompressedSeq's Links list.
func (cseq *CompressedSeq) Add(link *LinkToReference) {
	cseq.Links = append(cseq.Links, link)
}
