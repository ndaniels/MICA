package cablastp

import (
	"encoding/binary"
	"fmt"
	"os"
	"sync"
)

// Hard-coded file names for different pieces of a cablastp database.
const (
	FileCoarseFasta      = "coarse.fasta"
	FileCoarseFastaIndex = "coarse.fasta.index"
	FileCoarseLinks      = "coarse.links"
	FileCoarsePlainLinks = "coarse.links.plain"
	FileCoarseLinksIndex = "coarse.links.index"
	FileCoarseSeeds      = "coarse.seeds"
	FileCoarsePlainSeeds = "coarse.seeds.plain"
)

// CoarseDB represents a set of unique sequences that comprise the "coarse"
// database. Sequences in the coarse database, combined with information in the
// compressed database, are used to re-create the original sequences.
type CoarseDB struct {
	Seqs  []*CoarseSeq
	Seeds Seeds

	// The fastaCache is used during decompression. Namely, once a coarse
	// sequence is decompressed, it is cached into this map.
	fastaCache map[int]*CoarseSeq

	// The size of the coarse database index in bytes. This can be used to
	// quickly compute the number of sequences in the coarse database.
	// (Since each sequence is represented by a 64-bit integer offset, simply
	// divide by 8.)
	fastaIndexSize int64

	// File pointers to each file in the "coarse" part of a cablastp database.
	FileFasta      *os.File
	FileFastaIndex *os.File
	FileSeeds      *os.File
	FileLinks      *os.File
	FileLinksIndex *os.File

	// Ensures that adding a sequence to the coarse database is atomic.
	seqLock *sync.RWMutex

	// If a database is *created* without the read only flag set, then we have
	// to save the seeds table.
	// (We may want to remove this feature.)
	readOnly bool

	// If read only is not set, then this is used to track how many sequences
	// were already in the coarse database. (So that when we add more, we only
	// write the new ones.)
	seqsRead int

	// plain is a debugging feature that writes the links and seeds table (when
	// not read only) in a human readable format, rather than the default binary
	// format.
	plain bool

	// File pointers to use when 'plain' is true.
	plainLinks *os.File
	plainSeeds *os.File
}

// newWriteCoarseDB sets up a new coarse database to be written to (or opens
// an existing one ready for writing when 'appnd' is set).
func newWriteCoarseDB(appnd bool, db *DB) (*CoarseDB, error) {
	var err error

	Vprintln("\tOpening coarse database...")

	coarsedb := &CoarseDB{
		Seqs:           make([]*CoarseSeq, 0, 10000000),
		seqsRead:       0,
		Seeds:          NewSeeds(db.MapSeedSize, db.SeedLowComplexity),
		FileFasta:      nil,
		FileFastaIndex: nil,
		fastaIndexSize: 0,
		FileSeeds:      nil,
		FileLinks:      nil,
		FileLinksIndex: nil,
		seqLock:        &sync.RWMutex{},
		readOnly:       db.ReadOnly,
		plain:          db.SavePlain,
		plainSeeds:     nil,
	}
	coarsedb.FileFasta, err = db.openWriteFile(appnd, FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	coarsedb.FileFastaIndex, err = db.openWriteFile(appnd, FileCoarseFastaIndex)
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
	coarsedb.FileLinksIndex, err = db.openWriteFile(appnd, FileCoarseLinksIndex)
	if err != nil {
		return nil, err
	}

	info, err := coarsedb.FileFastaIndex.Stat()
	if err != nil {
		return nil, err
	}
	coarsedb.fastaIndexSize = info.Size()

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

		// After we've loaded the coarse database, the file offset should be
		// at the end of each file. For the coarse fasta file, this is
		// exactly what we want. But for the links and seeds files, we need
		// to clear the file and start over (since they are not amenable to
		// appending like the coarse fasta file is).
		// Do the same for plain files.
		trunc := func(f *os.File) (err error) {
			if err = f.Truncate(0); err != nil {
				return
			}
			if _, err = f.Seek(0, os.SEEK_SET); err != nil {
				return
			}
			return nil
		}
		if err = trunc(coarsedb.FileSeeds); err != nil {
			return nil, err
		}
		if err = trunc(coarsedb.FileLinks); err != nil {
			return nil, err
		}
		if err = trunc(coarsedb.FileLinksIndex); err != nil {
			return nil, err
		}
		if coarsedb.plain {
			if err = trunc(coarsedb.plainSeeds); err != nil {
				return nil, err
			}
			if err = trunc(coarsedb.plainLinks); err != nil {
				return nil, err
			}
		}
	}

	Vprintln("\tDone opening coarse database.")
	return coarsedb, nil
}

// newReadCoarseDB opens a coarse database and prepares it for reading. This
// is typically called before decompression.
func newReadCoarseDB(db *DB) (*CoarseDB, error) {
	var err error

	Vprintln("\tOpening coarse database...")

	coarsedb := &CoarseDB{
		Seqs:           make([]*CoarseSeq, 0, 100000),
		Seeds:          NewSeeds(db.MapSeedSize, db.SeedLowComplexity),
		FileFasta:      nil,
		fastaCache:     make(map[int]*CoarseSeq, 200),
		FileFastaIndex: nil,
		fastaIndexSize: 0,
		FileSeeds:      nil,
		FileLinks:      nil,
		FileLinksIndex: nil,
		seqLock:        nil,
		readOnly:       false,
		plain:          db.SavePlain,
	}
	coarsedb.FileFasta, err = db.openReadFile(FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	coarsedb.FileFastaIndex, err = db.openReadFile(FileCoarseFastaIndex)
	if err != nil {
		return nil, err
	}
	coarsedb.FileLinks, err = db.openReadFile(FileCoarseLinks)
	if err != nil {
		return nil, err
	}
	coarsedb.FileLinksIndex, err = db.openReadFile(FileCoarseLinksIndex)
	if err != nil {
		return nil, err
	}

	info, err := coarsedb.FileFastaIndex.Stat()
	if err != nil {
		return nil, err
	}
	coarsedb.fastaIndexSize = info.Size()

	Vprintln("\tDone opening coarse database.")
	return coarsedb, nil
}

// Add takes an original sequence, converts it to a coarse sequence, and
// adds it as a new coarse sequence to the coarse database. Seeds are
// also generated for each K-mer in the sequence. The resulting coarse
// sequence is returned along with its sequence identifier.
func (coarsedb *CoarseDB) Add(oseq []byte) (int, *CoarseSeq) {
	coarsedb.seqLock.Lock()
	id := len(coarsedb.Seqs)
	corSeq := NewCoarseSeq(id, "", oseq)
	coarsedb.Seqs = append(coarsedb.Seqs, corSeq)
	coarsedb.seqLock.Unlock()

	coarsedb.Seeds.Add(id, corSeq)

	return id, corSeq
}

// CoarseSeqGet is a thread-safe way to retrieve a sequence with index `i`
// from the coarse database.
func (coarsedb *CoarseDB) CoarseSeqGet(i uint) *CoarseSeq {
	coarsedb.seqLock.RLock()
	seq := coarsedb.Seqs[i]
	coarsedb.seqLock.RUnlock()

	return seq
}

// Expand will follow all links to compressed sequences for the coarse
// sequence at index `id` and return a slice of decompressed sequences.
func (coarsedb *CoarseDB) Expand(
	comdb *CompressedDB, id, start, end int) ([]OriginalSeq, error) {

	// Calculate the byte offset into the coarse links file where the links
	// for the coarse sequence `i` starts.
	off, err := coarsedb.linkOffset(id)
	if err != nil {
		return nil, fmt.Errorf("Could not get link offset: %s", err)
	}

	// Actually seek to that offset.
	newOff, err := coarsedb.FileLinks.Seek(off, os.SEEK_SET)
	if err != nil {
		return nil, fmt.Errorf("Could not seek: %s", err)
	} else if newOff != off {
		return nil,
			fmt.Errorf("Tried to seek to offset %d in the coarse links, "+
				"but seeked to %d instead.", off, newOff)
	}

	// Read in the number of links for this sequence.
	// Each link corresponds to a single original sequence.
	var numLinks uint32
	err = binary.Read(coarsedb.FileLinks, binary.BigEndian, &numLinks)
	if err != nil {
		return nil, fmt.Errorf("Could not read number of links: %s", err)
	}

	// We use a map as a set of original sequence ids for eliminating
	// duplicates (since a coarse sequence can point to different pieces of the
	// same compressed sequence).
	ids := make(map[uint32]bool, numLinks)
	oseqs := make([]OriginalSeq, 0, numLinks)
	s, e := uint16(start), uint16(end)
	for i := uint32(0); i < numLinks; i++ {
		compLink, err := coarsedb.readLink()
		if err != nil {
			return nil, fmt.Errorf("Could not read link: %s", err)
		}

		// We only use this link if the match is in the range.
		if e < compLink.CoarseStart || s > compLink.CoarseEnd {
			continue
		}

		// Don't decompress the same original sequence more than once.
		if ids[compLink.OrgSeqId] {
			continue
		}

		oseq, err := comdb.ReadSeq(coarsedb, int(compLink.OrgSeqId))
		if err != nil {
			return nil, fmt.Errorf(
				"Could not read compressed sequence: %s", err)
		}
		ids[compLink.OrgSeqId] = true
		oseqs = append(oseqs, oseq)
	}

	return oseqs, nil
}

// NumRequences returns the number of sequences in the coarse database based
// on the file size of the coarse database index.
func (coarsedb *CoarseDB) NumSequences() int {
	return int(coarsedb.fastaIndexSize / 8)
}

// ReadCoarseSeq reads the coarse sequence with identifier 'id' from disk, using
// the fasta index. (If a coarse sequence has already been read, it is returned
// from cache to save trips to disk.)
//
// TODO: Note that this does *not* recover links typically found in a coarse
// sequence, although it probably should to avoid doing it in CoarseDB.Expand.
func (coarsedb *CoarseDB) ReadCoarseSeq(id int) (*CoarseSeq, error) {
	// Prevent reading the same coarse sequence over and over.
	if coarseSeq, ok := coarsedb.fastaCache[id]; ok {
		return coarseSeq, nil
	}

	off, err := coarsedb.coarseOffset(id)
	if err != nil {
		return nil, fmt.Errorf("Could not get coarse offset: %s", err)
	}

	newOff, err := coarsedb.FileFasta.Seek(off, os.SEEK_SET)
	if err != nil {
		return nil, fmt.Errorf("Could not seek in coarse fasta: %s", err)
	} else if newOff != off {
		return nil,
			fmt.Errorf("Tried to seek to offset %d in the coarse fasta file, "+
				"but seeked to %d instead.", off, newOff)
	}

	// Read in the sequence.
	var corSeqId int
	var residues string
	n, err := fmt.Fscanf(coarsedb.FileFasta, "> %d\n%s\n", &corSeqId, &residues)
	if err != nil {
		return nil, fmt.Errorf("Could not scan coarse sequence %d: %s", id, err)
	} else if n != 2 {
		return nil, fmt.Errorf("Expected to read in two values for coarse "+
			"sequence %d, but read %d values instead.", id, n)
	} else if corSeqId != id {
		return nil, fmt.Errorf("Expected to read coarse sequence %d but read "+
			"coarse sequence %d instead.", id, corSeqId)
	}

	coarseSeq := NewCoarseSeq(id, "", []byte(residues))
	coarsedb.fastaCache[id] = coarseSeq
	return coarseSeq, nil
}

// coarseOffset returns the integer byte offset into the coarse database of
// a particular coarse sequence. The offset is read from the coarse database
// index.
//
// An error is returned if the file seek fails.
func (coarsedb *CoarseDB) coarseOffset(id int) (seqOff int64, err error) {
	tryOff := int64(id) * 8
	realOff, err := coarsedb.FileFastaIndex.Seek(tryOff, os.SEEK_SET)
	if err != nil {
		return
	} else if tryOff != realOff {
		return 0,
			fmt.Errorf("Tried to seek to offset %d in the coarse index, "+
				"but seeked to %d instead.", tryOff, realOff)
	}
	err = binary.Read(coarsedb.FileFastaIndex, binary.BigEndian, &seqOff)
	return
}

// linkOffset returns the integer byte offset into the coarse links database
// of a particular coarse sequence. The offset is read from the coarse links
// database index.
//
// An error is returned if the file seek fails.
func (coarsedb *CoarseDB) linkOffset(id int) (seqOff int64, err error) {
	tryOff := int64(id) * 8
	realOff, err := coarsedb.FileLinksIndex.Seek(tryOff, os.SEEK_SET)
	if err != nil {
		return
	} else if tryOff != realOff {
		return 0,
			fmt.Errorf("Tried to seek to offset %d in the coarse links index, "+
				"but seeked to %d instead.", tryOff, realOff)
	}
	err = binary.Read(coarsedb.FileLinksIndex, binary.BigEndian, &seqOff)
	return
}

// readClose closes all files necessary for reading the coarse database.
func (coarsedb *CoarseDB) readClose() {
	coarsedb.FileFasta.Close()
	coarsedb.FileFastaIndex.Close()
	coarsedb.FileLinks.Close()
	coarsedb.FileLinksIndex.Close()
}

// writeClose closes all files necessary for writing the coarse database.
func (coarsedb *CoarseDB) writeClose() {
	coarsedb.FileFasta.Close()
	coarsedb.FileFastaIndex.Close()
	coarsedb.FileSeeds.Close()
	coarsedb.FileLinks.Close()
	coarsedb.FileLinksIndex.Close()
	if coarsedb.plain {
		coarsedb.plainLinks.Close()
		coarsedb.plainSeeds.Close()
	}
}

// load reads the entire coarse database (sequences and links) into memory.
// If the database is being appended to, the seeds table is also read into
// memory.
//
// (This is only called when a coarse database is being appended to.)
func (coarsedb *CoarseDB) load() (err error) {
	if err = coarsedb.readFasta(); err != nil {
		return
	}
	if err = coarsedb.readLinks(); err != nil {
		return
	}
	if coarsedb.FileSeeds != nil {
		if err = coarsedb.readSeeds(); err != nil {
			return
		}
	}
	return nil
}

// save will save the coarse database as a FASTA file and a binary
// encoding of all coarse links. Seeds are also saved if this is not a read
// only database.
func (coarsedb *CoarseDB) save() error {
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
