package cablastp

import (
	"bytes"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"os"
)

type CompressedDB struct {
	SeqCache   map[int]CompressedSeq
	File       *os.File
	Index      *os.File
	writerChan chan CompressedSeq
	writerDone chan struct{}
}

func NewCompressedDB(file *os.File, index *os.File) *CompressedDB {
	cdb := &CompressedDB{
		SeqCache:   nil,
		File:       file,
		Index:      index,
		writerChan: make(chan CompressedSeq, 500),
		writerDone: make(chan struct{}, 0),
	}
	go cdb.writer()

	return cdb
}

func NewAppendCompressedDB(file *os.File, index *os.File) *CompressedDB {
	return NewCompressedDB(file, index)
}

func LoadCompressedDB(file *os.File, index *os.File) (*CompressedDB, error) {
	cdb := &CompressedDB{
		SeqCache:   make(map[int]CompressedSeq, 100),
		File:       file,
		Index:      index,
		writerChan: nil,
		writerDone: nil,
	}
	return cdb, nil
}

func (comdb *CompressedDB) ReadClose() {
	comdb.File.Close()
	comdb.Index.Close()
}

func (comdb *CompressedDB) WriteClose() {
	close(comdb.writerChan) // will close comdb.File

	// Wait for the writer goroutine to finish.
	<-comdb.writerDone

	comdb.Index.Close()
}

func (comdb *CompressedDB) LoadSeq(orgSeqId int) *CompressedSeq {
	if comdb.writerChan != nil {
		panic(fmt.Sprintf("A compressed database cannot be read while it is " +
			"also being modified."))
	}
	return nil
}

func (comdb *CompressedDB) Write(cseq CompressedSeq) {
	comdb.writerChan <- cseq
}

func (comdb *CompressedDB) writer() {
	var record []string
	var err error

	byteOffset := int64(0)
	buf := new(bytes.Buffer)
	csvWriter := csv.NewWriter(buf)
	csvWriter.Comma = ','
	csvWriter.UseCRLF = false

	for cseq := range comdb.writerChan {
		// Reset the buffer so it's empty. We want it to only contain
		// the next record we're writing.
		buf.Reset()

		// Allocate memory for creating the next record.
		// A record is a sequence name followed by four-tuples of links:
		// (coarse-seq-id, coarse-start, coarse-end, diff).
		record = make([]string, 0, 1+4*len(cseq.Links))
		record = append(record, cseq.Name)
		for _, link := range cseq.Links {
			record = append(record,
				fmt.Sprintf("%d", link.CoarseSeqId),
				fmt.Sprintf("%d", link.CoarseStart),
				fmt.Sprintf("%d", link.CoarseEnd),
				link.Diff)
		}

		// Write the record to our *buffer* and flush it.
		if err = csvWriter.Write(record); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}
		csvWriter.Flush()

		// Pass the bytes on to the compressed file.
		if _, err = comdb.File.Write(buf.Bytes()); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}

		// Now write the byte offset that points to the start of this record.
		err = binary.Write(comdb.Index, binary.BigEndian, byteOffset)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}

		// Increment the byte offset to be at the end of this record.
		byteOffset += int64(buf.Len())
	}
	comdb.File.Close()
	comdb.writerDone <- struct{}{}
}

// CompressedSeq corresponds to the components of a compressed sequence.
type CompressedSeq struct {
	Id int

	// Name is an uncompressed string from the original FASTA header.
	Name string

	// Links is an ordered lists of links to portions of the reference
	// database. When all links are followed, the concatenation of each
	// sequence correspond to each link equals the entire original sequence.
	Links []LinkToCoarse
}

// NewCompressedSeq creates a CompressedSeq value using the name provided.
// The Link slice is initialized but empty.
func NewCompressedSeq(id int, name string) CompressedSeq {
	return CompressedSeq{
		Id:    id,
		Name:  name,
		Links: make([]LinkToCoarse, 0, 10),
	}
}

// Add will add a LinkToReference to the end of the CompressedSeq's Links list.
func (cseq *CompressedSeq) Add(link LinkToCoarse) {
	cseq.Links = append(cseq.Links, link)
}
