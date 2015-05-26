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

// A CompressedDB corresponds to a list of all original sequences compressed
// by replacing regions of sequences that are redundant with pointers to
// similar regions in the coarse database. Each pointer includes an offset and
// an edit script, which allows complete recovery of the original sequence.
//
// N.B. A compressed database doesn't keep an in memory representation of
// all compressed sequences. In particular, writing to a compressed database
// always corresponds to writing a compressed sequence to disk. And reading
// from a compressed database always corresponds to reading a sequence from
// disk (unless it has been cached in 'seqCache').
type CompressedDB struct {
	// File pointers to be used in reading/writing compressed databases.
	File  *os.File
	Index *os.File

	// The size of the compressed database index in bytes. Since the index
	// contains precisely one 64-bit integer byte offset for every sequence
	// in the compressed database, the index size can be used to quickly
	// compute the number of sequences in the compressed database.
	indexSize int64

	// A pair of channels used to facilitate writing compressed sequences as
	// they are processed during compression. (The writer operates in its own
	// gorotuine.)
	writerChan chan CompressedSeq
	writerDone chan struct{}

	// A compressed database is stored in CSV format. Each CSV record contains
	// the original sequence's header, followed by a list of quadruples, where
	// each quadruple is a pointer to a region in the coarse database: a coarse
	// sequence identifier, the start/end of the coarse sequence, and an edit
	// script. Combined, this information can recover the original sequence
	// in full.
	csvReader *csv.Reader

	// Caches already read sequences from the compressed database while reading.
	seqCache map[int]OriginalSeq
}

// newWriteCompressedDB creates a new compressed database ready for writing
// (or opens an existing compressed database and prepares it for appending).
func newWriteCompressedDB(appnd bool, db *DB) (*CompressedDB, error) {
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

// newReadCompressedDB opens a compressed database and prepares it for reading.
func newReadCompressedDB(db *DB) (*CompressedDB, error) {
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

	info, err := cdb.Index.Stat()
	if err != nil {
		return nil, err
	}
	cdb.indexSize = info.Size()

	cdb.csvReader = csv.NewReader(cdb.File)
	cdb.csvReader.Comma = ','
	cdb.csvReader.FieldsPerRecord = -1

	Vprintln("\tDone opening compressed database.")
	return cdb, nil
}

// SeqGet reads a sequence from the compressed database, and decompressed it
// using the coarse database provided. The decompressed sequence is then added
// to cache.
//
// If the sequence has already been decompressed, the decompressed sequence
// from cache is returned.
//
// SeqGet will panic if it is called while a compressed database is open for
// writing.
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

// NumSequences returns the number of sequences in the compressed database
// using the file size of the index.
func (comdb *CompressedDB) NumSequences() int {
	return int(comdb.indexSize / 8)
}

// readClose closes all appropriate files used in reading a compressed database.
func (comdb *CompressedDB) readClose() {
	comdb.File.Close()
	comdb.Index.Close()
}

// writeClose closes all appropriate files used in writing a compressed
// database. It will also wait for the writer goroutine to finish.
func (comdb *CompressedDB) writeClose() {
	close(comdb.writerChan) // will close comdb.File

	// Wait for the writer goroutine to finish.
	<-comdb.writerDone
}

// Write queues a new compressed sequence to be written to disk.
func (comdb *CompressedDB) Write(cseq CompressedSeq) {
	comdb.writerChan <- cseq
}

// CompressedSeq corresponds to the components of a compressed sequence.
type CompressedSeq struct {
	// A sequence number.
	Id int

	// Name is an uncompressed string from the original FASTA header.
	Name string

	// Links is an ordered lists of links to portions of the reference
	// database. When all links are followed, the concatenation of each
	// sequence corresponding to each link equals the entire original sequence.
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

// Add will add a LinkToCoarse to the end of the CompressedSeq's Links list.
func (cseq *CompressedSeq) Add(link LinkToCoarse) {
	cseq.Links = append(cseq.Links, link)
}

// Decompress decompresses a particular compressed sequence using the given
// coarse sequence. Namely, all of the links are followed and all of the
// edit scripts are "applied" to recover the original sequence.
func (cseq CompressedSeq) Decompress(coarse *CoarseDB) (OriginalSeq, error) {
	residues := make([]byte, 0, 20)
	for _, lk := range cseq.Links {
		if lk.CoarseSeqId < 0 || lk.CoarseSeqId >= uint(coarse.NumSequences()) {
			return OriginalSeq{},
				fmt.Errorf("Cannot decompress compressed sequence (id: %d), "+
					"because a link refers to an invalid coarse sequence "+
					"id: %d.", cseq.Id, lk.CoarseSeqId)
		}
		editScript, err := NewEditScriptParse(lk.Diff)
		if err != nil {
			return OriginalSeq{}, err
		}

		coarseSeq, err := coarse.ReadCoarseSeq(int(lk.CoarseSeqId))
		if err != nil {
			return OriginalSeq{}, err
		}
		subCorres := coarseSeq.Residues[lk.CoarseStart:lk.CoarseEnd]
		residues = append(residues, editScript.Apply(subCorres)...)
	}
	return *NewOriginalSeq(cseq.Id, cseq.Name, residues), nil
}
