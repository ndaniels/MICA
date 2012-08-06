package cablastp

import (
	"bufio"
	"compress/gzip"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"os"
	"sync"
)

const (
	FileCoarseFasta      = "coarse.fasta"
	FileCoarseLinks      = "coarse.links"
	FileCoarsePlainLinks = "coarse.links.plain"
	FileCoarseSeeds      = "coarse.seeds"
	FileCoarsePlainSeeds = "coarse.seeds.plain"
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

	plain      bool
	plainLinks *os.File
	plainSeeds *os.File
}

// NewCoarseDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func NewWriteCoarseDB(appnd bool, db *DB) (*CoarseDB, error) {
	var err error

	coarsedb := &CoarseDB{
		Seqs:       make([]*CoarseSeq, 0, 10000000),
		Seeds:      NewSeeds(db.MapSeedSize),
		FileFasta:  nil,
		FileSeeds:  nil,
		FileLinks:  nil,
		seqLock:    &sync.RWMutex{},
		plain:      db.SavePlain,
		plainSeeds: nil,
	}
	coarsedb.FileFasta, err = db.openWriteFile(false, FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	coarsedb.FileSeeds, err = db.openWriteFile(false, FileCoarseSeeds)
	if err != nil {
		return nil, err
	}
	coarsedb.FileLinks, err = db.openWriteFile(false, FileCoarseLinks)
	if err != nil {
		return nil, err
	}

	if coarsedb.plain {
		coarsedb.plainLinks, err = db.openWriteFile(false, FileCoarsePlainLinks)
		if err != nil {
			return nil, err
		}
		coarsedb.plainSeeds, err = db.openWriteFile(false, FileCoarsePlainSeeds)
		if err != nil {
			return nil, err
		}
	}
	return coarsedb, nil
}

func NewReadCoarseDB(db *DB) (*CoarseDB, error) {
	var err error

	coarsedb := &CoarseDB{
		Seqs:      make([]*CoarseSeq, 0, 10000000),
		Seeds:     NewSeeds(db.MapSeedSize),
		FileFasta: nil,
		FileSeeds: nil,
		FileLinks: nil,
		seqLock:   nil,
		plain:     db.SavePlain,
	}
	coarsedb.FileFasta, err = db.openReadFile(FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	coarsedb.FileLinks, err = db.openReadFile(FileCoarseLinks)
	if err != nil {
		return nil, err
	}
	return coarsedb, nil
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
	if coarsedb.plain {
		coarsedb.plainLinks.Close()
		coarsedb.plainSeeds.Close()
	}
}

// Save will save the reference database as a coarse FASTA file and a binary
// encoding of all reference links.
func (coarsedb *CoarseDB) Save() error {
	coarsedb.seqLock.RLock()
	defer coarsedb.seqLock.RUnlock()

	if err := coarsedb.saveFasta(); err != nil {
		return err
	}
	if err := coarsedb.saveLinks(); err != nil {
		return err
	}
	if err := coarsedb.saveSeeds(); err != nil {
		return err
	}
	if coarsedb.plain {
		if err := coarsedb.saveLinksPlain(); err != nil {
			return err
		}
		if err := coarsedb.saveSeedsPlain(); err != nil {
			return err
		}
	}
	return nil
}

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

func (coarsedb *CoarseDB) saveLinks() error {
	gzipWriter, _ := gzip.NewWriterLevel(coarsedb.FileLinks, gzip.BestSpeed)
	for _, seq := range coarsedb.Seqs {
		for link := seq.Links; link != nil; link = link.Next {
			binary.Write(gzipWriter, binary.BigEndian, link.OrgSeqId)
			binary.Write(gzipWriter, binary.BigEndian, link.CoarseStart)
			binary.Write(gzipWriter, binary.BigEndian, link.CoarseEnd)
		}
		if _, err := gzipWriter.Write([]byte{'\n'}); err != nil {
			return err
		}
	}
	if err := gzipWriter.Close(); err != nil {
		return nil
	}
	return nil
}

func (coarsedb *CoarseDB) saveLinksPlain() error {
	csvWriter := csv.NewWriter(coarsedb.plainLinks)
	record := make([]string, 0, 10)
	for _, seq := range coarsedb.Seqs {
		record = record[:0]
		for link := seq.Links; link != nil; link = link.Next {
			record = append(record,
				fmt.Sprintf("%d", link.OrgSeqId),
				fmt.Sprintf("%d", link.CoarseStart),
				fmt.Sprintf("%d", link.CoarseEnd))
		}
		if err := csvWriter.Write(record); err != nil {
			return err
		}
	}
	csvWriter.Flush()
	return nil
}

func (coarsedb *CoarseDB) saveSeeds() error {
	var i int32

	gzipWriter, _ := gzip.NewWriterLevel(coarsedb.FileSeeds, gzip.BestSpeed)
	for i = 0; i < int32(coarsedb.Seeds.powers[coarsedb.Seeds.SeedSize]); i++ {
		if coarsedb.Seeds.Locs[i] == nil {
			continue
		}

		binary.Write(gzipWriter, binary.BigEndian, i)
		for loc := coarsedb.Seeds.Locs[i]; loc != nil; loc = loc.Next {
			binary.Write(gzipWriter, binary.BigEndian, loc.SeqInd)
			binary.Write(gzipWriter, binary.BigEndian, loc.ResInd)
		}
		if _, err := gzipWriter.Write([]byte{'\n'}); err != nil {
			return err
		}
	}
	if err := gzipWriter.Close(); err != nil {
		return err
	}
	return nil
}

func (coarsedb *CoarseDB) saveSeedsPlain() error {
	csvWriter := csv.NewWriter(coarsedb.plainSeeds)
	record := make([]string, 0, 10)
	for i := 0; i < coarsedb.Seeds.powers[coarsedb.Seeds.SeedSize]; i++ {
		if coarsedb.Seeds.Locs[i] == nil {
			continue
		}

		record = record[:0]
		record = append(record, string(coarsedb.Seeds.unhashKmer(i)))
		for loc := coarsedb.Seeds.Locs[i]; loc != nil; loc = loc.Next {
			record = append(record,
				fmt.Sprintf("%d", loc.SeqInd),
				fmt.Sprintf("%d", loc.ResInd))
		}
		if err := csvWriter.Write(record); err != nil {
			return err
		}
	}
	csvWriter.Flush()
	return nil
}
