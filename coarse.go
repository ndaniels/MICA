package cablastp

import (
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
	Seqs     []*CoarseSeq
	seqsRead int
	Seeds    Seeds

	FileFasta *os.File
	FileSeeds *os.File
	FileLinks *os.File

	seqLock *sync.RWMutex

	readOnly   bool
	plain      bool
	plainLinks *os.File
	plainSeeds *os.File
}

// NewCoarseDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func NewWriteCoarseDB(appnd bool, db *DB) (*CoarseDB, error) {
	var err error

	Vprintln("\tOpening coarse database...")

	coarsedb := &CoarseDB{
		Seqs:       make([]*CoarseSeq, 0, 10000000),
		seqsRead:   0,
		Seeds:      NewSeeds(db.MapSeedSize),
		FileFasta:  nil,
		FileSeeds:  nil,
		FileLinks:  nil,
		seqLock:    &sync.RWMutex{},
		readOnly:   db.ReadOnly,
		plain:      db.SavePlain,
		plainSeeds: nil,
	}
	coarsedb.FileFasta, err = db.openWriteFile(appnd, FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	coarsedb.FileSeeds, err = db.openWriteFile(appnd, FileCoarseSeeds)
	if err != nil {
		return nil, err
	}
	coarsedb.FileLinks, err = db.openWriteFile(appnd, FileCoarseLinks)
	if err != nil {
		return nil, err
	}

	if coarsedb.plain {
		coarsedb.plainLinks, err = db.openWriteFile(appnd, FileCoarsePlainLinks)
		if err != nil {
			return nil, err
		}
		coarsedb.plainSeeds, err = db.openWriteFile(appnd, FileCoarsePlainSeeds)
		if err != nil {
			return nil, err
		}
	}

	if appnd {
		if err = coarsedb.load(); err != nil {
			return nil, err
		}
	}

	Vprintln("\tDone opening coarse database.")
	return coarsedb, nil
}

func NewReadCoarseDB(db *DB) (*CoarseDB, error) {
	var err error

	Vprintln("\tOpening coarse database...")

	coarsedb := &CoarseDB{
		Seqs:      make([]*CoarseSeq, 0, 10000000),
		Seeds:     NewSeeds(db.MapSeedSize),
		FileFasta: nil,
		FileSeeds: nil,
		FileLinks: nil,
		seqLock:   nil,
		readOnly:  false,
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

	if err := coarsedb.load(); err != nil {
		return nil, err
	}

	Vprintln("\tDone opening coarse database.")
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

func (coarsedb *CoarseDB) load() (err error) {
	if err = coarsedb.readFasta(); err != nil {
		return
	}
	if err = coarsedb.readSeeds(); err != nil {
		return
	}
	if err = coarsedb.readLinks(); err != nil {
		return
	}

	// After we've loaded the coarse database, the file offset should be
	// at the end of each file. For the coarse fasta file, this is
	// exactly what we want. But for the links and seeds files, we need
	// to clear the file and start over (since it is not amenable to
	// appending like the coarse fasta file is).
	// Do the same for plain files.
	if err = coarsedb.FileSeeds.Truncate(0); err != nil {
		return
	}
	if _, err = coarsedb.FileSeeds.Seek(0, os.SEEK_SET); err != nil {
		return
	}

	if err = coarsedb.FileLinks.Truncate(0); err != nil {
		return
	}
	if _, err = coarsedb.FileLinks.Seek(0, os.SEEK_SET); err != nil {
		return
	}

	if coarsedb.plain {
		if err = coarsedb.plainSeeds.Truncate(0); err != nil {
			return
		}
		if _, err = coarsedb.plainSeeds.Seek(0, os.SEEK_SET); err != nil {
			return
		}

		if err = coarsedb.plainLinks.Truncate(0); err != nil {
			return
		}
		if _, err = coarsedb.plainLinks.Seek(0, os.SEEK_SET); err != nil {
			return
		}
	}

	return nil
}

// Save will save the reference database as a coarse FASTA file and a binary
// encoding of all reference links.
func (coarsedb *CoarseDB) Save() error {
	coarsedb.seqLock.RLock()
	defer coarsedb.seqLock.RUnlock()

	errc := make(chan error, 20)
	wg := &sync.WaitGroup{}

	wg.Add(1)
	go func() {
		if err := coarsedb.saveFasta(); err != nil {
			errc <- err
		}
		wg.Done()
	}()

	wg.Add(1)
	go func() {
		if err := coarsedb.saveLinks(); err != nil {
			errc <- err
		}
		wg.Done()
	}()

	if !coarsedb.readOnly {
		wg.Add(1)
		go func() {
			if err := coarsedb.saveSeeds(); err != nil {
				errc <- err
			}
			wg.Done()
		}()
	}
	if coarsedb.plain {
		wg.Add(1)
		go func() {
			if err := coarsedb.saveLinksPlain(); err != nil {
				errc <- err
			}
			wg.Done()
		}()
		if !coarsedb.readOnly {
			wg.Add(1)
			go func() {
				if err := coarsedb.saveSeedsPlain(); err != nil {
					errc <- err
				}
				wg.Done()
			}()
		}
	}
	wg.Wait()

	// If there's something in the error channel, pop off the first
	// error and return that.
	if len(errc) > 0 {
		return <-errc
	}
	return nil
}
