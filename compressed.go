package cablastp

import (
	"encoding/csv"
	"fmt"
	"os"
	"strings"
)

const (
	FileCompressed = "compressed"
	FileIndex      = "compressed.index"
)

type CompressedDB struct {
	seqCache   map[int]OriginalSeq
	File       *os.File
	Index      *os.File
	indexSize  int64
	writerChan chan CompressedSeq
	writerDone chan struct{}
	csvReader  *csv.Reader
}

func NewWriteCompressedDB(appnd bool, db *DB) (*CompressedDB, error) {
	var err error

	Vprintln("\tOpening compressed database...")

	cdb := &CompressedDB{
		seqCache:   nil,
		File:       nil,
		Index:      nil,
		writerChan: make(chan CompressedSeq, 500),
		writerDone: make(chan struct{}, 0),
	}

	fileFlags := os.O_RDWR | os.O_CREATE | os.O_TRUNC
	if appnd {
		fileFlags = os.O_RDWR | os.O_APPEND
	}
	cdb.File, err = os.OpenFile(db.filePath(FileCompressed), fileFlags, 0666)
	if err != nil {
		return nil, err
	}
	cdb.Index, err = os.OpenFile(db.filePath(FileIndex), fileFlags, 0666)
	if err != nil {
		return nil, err
	}

	info, err := cdb.Index.Stat()
	if err != nil {
		return nil, err
	}
	cdb.indexSize = info.Size()

	go cdb.writer()

	Vprintln("\tDone opening compressed database.")
	return cdb, nil
}

func NewReadCompressedDB(db *DB) (*CompressedDB, error) {
	var err error

	Vprintln("\tOpening compressed database...")

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

	Vprintln("\tDone opening compressed database.")
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
		comdb.seqCache[orgSeqId], err = comdb.ReadSeq(coarsedb, orgSeqId)
		if err != nil {
			return OriginalSeq{}, err
		}
	}
	return comdb.seqCache[orgSeqId], nil
}

func (comdb *CompressedDB) NumSequences() int {
	return int(comdb.indexSize / 8)
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

func (cseq CompressedSeq) String() string {
	lines := make([]string, len(cseq.Links))
	for i, link := range cseq.Links {
		lines[i] = fmt.Sprintf("coarse id: %d, start: %d, end: %d\n%s",
			link.CoarseSeqId, link.CoarseStart, link.CoarseEnd, link.Diff)
	}
	return strings.Join(lines, "\n")
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
