package cablastp

import (
	"bytes"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
)

const (
	FileCompressed = "compressed.cbp"
	FileIndex      = "index"
)

type CompressedDB struct {
	seqCache   map[int]OriginalSeq
	File       *os.File
	Index      *os.File
	writerChan chan CompressedSeq
	writerDone chan struct{}
	csvReader  *csv.Reader
}

func NewWriteCompressedDB(appnd bool, db *DB) (*CompressedDB, error) {
	var err error

	cdb := &CompressedDB{
		seqCache:   nil,
		File:       nil,
		Index:      nil,
		writerChan: make(chan CompressedSeq, 500),
		writerDone: make(chan struct{}, 0),
	}
	cdb.File, err = db.openWriteFile(appnd, FileCompressed)
	if err != nil {
		return nil, err
	}
	cdb.Index, err = db.openWriteFile(appnd, FileIndex)
	if err != nil {
		return nil, err
	}

	go cdb.writer()

	return cdb, nil
}

func NewReadCompressedDB(db *DB) (*CompressedDB, error) {
	var err error

	cdb := &CompressedDB{
		seqCache:   make(map[int]OriginalSeq, 100),
		File:       nil,
		Index:      nil,
		writerChan: nil,
		writerDone: nil,
		csvReader:  nil,
	}
	cdb.File, err = db.openReadFile(FileCompressed)
	if err != nil {
		return nil, err
	}
	cdb.Index, err = db.openReadFile(FileIndex)
	if err != nil {
		return nil, err
	}
	cdb.csvReader = csv.NewReader(cdb.File)
	return cdb, nil
}

func (comdb *CompressedDB) ReadClose() {
	comdb.File.Close()
	comdb.Index.Close()
}

func (comdb *CompressedDB) SeqGet(
	coarsedb *CoarseDB, orgSeqId int) (OriginalSeq, error) {

	var err error

	if comdb.writerChan != nil {
		panic(fmt.Sprintf("A compressed database cannot be read while it is " +
			"also being modified."))
	}
	if _, ok := comdb.seqCache[orgSeqId]; !ok {
		comdb.seqCache[orgSeqId], err = comdb.readSeq(coarsedb, orgSeqId)
		if err != nil {
			return OriginalSeq{}, err
		}
	}
	return comdb.seqCache[orgSeqId], nil
}

func (comdb *CompressedDB) readSeq(
	coarsedb *CoarseDB, orgSeqId int) (OriginalSeq, error) {

	off, err := comdb.orgSeqOffset(orgSeqId)
	if err != nil {
		return OriginalSeq{}, err
	}

	newOff, err := comdb.File.Seek(off, 0)
	if err != nil {
		return OriginalSeq{}, err
	} else if newOff != off {
		return OriginalSeq{},
			fmt.Errorf("Tried to seek to offset %d in the compressed "+
				"database, but seeked to %d instead.", off, newOff)
	}

	record, err := comdb.csvReader.Read()
	if err != nil {
		return OriginalSeq{}, err
	}

	cseq, err := readCompressedSeq(orgSeqId, record)
	if err != nil {
		return OriginalSeq{}, err
	}
	return cseq.Decompress(coarsedb)
}

func (comdb *CompressedDB) orgSeqOffset(id int) (seqOff int64, err error) {
	tryOff := int64(id) * 8
	realOff, err := comdb.Index.Seek(tryOff, 0)
	if err != nil {
		return 0, err
	} else if tryOff != realOff {
		return 0,
			fmt.Errorf("Tried to seek to offset %d in the compressed index, "+
				"but seeked to %d instead.", tryOff, realOff)
	}
	err = binary.Read(comdb.Index, binary.BigEndian, &seqOff)
	return
}

func (comdb *CompressedDB) WriteClose() {
	close(comdb.writerChan) // will close comdb.File

	// Wait for the writer goroutine to finish.
	<-comdb.writerDone

	comdb.Index.Close()
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

func readCompressedSeq(id int, record []string) (CompressedSeq, error) {
	cseq := CompressedSeq{
		Id:    id,
		Name:  string([]byte(record[0])),
		Links: make([]LinkToCoarse, (len(record)-1)/4),
	}

	for i := 1; i < len(record); i += 4 {
		coarseSeqId64, err := strconv.Atoi(record[i+0])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseStart64, err := strconv.Atoi(record[i+1])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseEnd64, err := strconv.Atoi(record[i+2])
		if err != nil {
			return CompressedSeq{}, nil
		}
		lk := NewLinkToCoarseNoDiff(
			int(coarseSeqId64), int(coarseStart64), int(coarseEnd64))
		lk.Diff = string([]byte(record[i+3]))

		cseq.Add(lk)
	}
	return cseq, nil
}

// Add will add a LinkToReference to the end of the CompressedSeq's Links list.
func (cseq *CompressedSeq) Add(link LinkToCoarse) {
	cseq.Links = append(cseq.Links, link)
}

func (cseq CompressedSeq) Decompress(coarsedb *CoarseDB) (OriginalSeq, error) {
	var corres, subCorres []byte
	residues := make([]byte, 0, 20)
	for _, lk := range cseq.Links {
		if lk.CoarseSeqId < 0 || lk.CoarseSeqId >= len(coarsedb.Seqs) {
			return OriginalSeq{},
				fmt.Errorf("Cannot decompress compressed sequence (id: %d), "+
					"because a link refers to an invalid coarse sequence "+
					"id: %d.", cseq.Id, lk.CoarseSeqId)
		}
		editScript, err := NewEditScriptParse(lk.Diff)
		if err != nil {
			return OriginalSeq{}, err
		}

		corres = coarsedb.Seqs[lk.CoarseSeqId].Residues
		subCorres = corres[lk.CoarseStart:lk.CoarseEnd]
		residues = append(residues, editScript.Apply(subCorres)...)
	}
	return *NewOriginalSeq(cseq.Id, cseq.Name, residues), nil
}
