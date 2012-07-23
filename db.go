package cablastp

import (
	"fmt"
	"os"
	"path"
)

const (
	FileCoarseFasta = "coarse.fasta"
	FileCoarseLinks = "coarse.links"
	FileCoarseSeeds = "coarse.seeds"
	FileCompressed  = "compressed.cbp"
	FileIndex       = "index"
)

type DB struct {
	Name     string
	ComDB    *CompressedDB
	CoarseDB *CoarseDB

	coarseFasta, coarseSeeds, coarseLinks, compressed, index *os.File
}

func NewDB(dir string, seedSize int, add, plain bool) (*DB, error) {
	// If we're not adding to a database, check to make sure the directory
	// doesn't exist.
	if !add {
		_, err := os.Open(dir)
		if err == nil {
			return nil, fmt.Errorf("The directory '%s' already exists. A "+
				"new compressed database cannot be created in the same "+
				"directory as an existing database. If you want to append to "+
				"to an existing database with, use the '--append' flag.", dir)
		}
		if err != nil && !os.IsNotExist(err) {
			return nil, fmt.Errorf("An error occurred when checking if '%s' "+
				"exists: %s.", dir, err)
		}
	} else { // otherwise, check to make sure it *does* exist.
		_, err := os.Open(dir)
		if err != nil {
			return nil, fmt.Errorf("Could not open '%s' for appending "+
				"because: %s.", dir, err)
		}
	}

	// If we're not adding, make the directory.
	err := os.Mkdir(dir, 0777)
	if err != nil {
		return nil, fmt.Errorf("Could not create directory '%s': %s.", dir, err)
	}

	// Initialize the DB struct. We add the file pointers and databases
	// later.
	db := &DB{
		Name: path.Base(dir),
	}

	// Now we need to open five files:
	// 1) Coarse FASTA database
	// 2) Coarse links
	// 3) Seeds
	// 4) Original links
	// 5) A line index
	//
	// If 'add' is false, then each file is opened and truncated for writing.
	// If 'add' is true, then each file is opened in 'append' mode for writing.
	open := func(name string) (*os.File, error) {
		return openDbFile(add, dir, name)
	}

	db.coarseFasta, err = open(FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	db.coarseLinks, err = open(FileCoarseLinks)
	if err != nil {
		return nil, err
	}
	db.coarseSeeds, err = open(FileCoarseSeeds)
	if err != nil {
		return nil, err
	}
	db.compressed, err = open(FileCompressed)
	if err != nil {
		return nil, err
	}
	db.index, err = open(FileIndex)
	if err != nil {
		return nil, err
	}

	db.ComDB = NewCompressedDB(db.compressed, db.index, plain)
	db.CoarseDB = NewCoarseDB(db.coarseFasta, db.coarseSeeds, db.coarseLinks,
		seedSize, plain)

	return db, nil
}

func (db *DB) Close() {
	db.CoarseDB.Close()
	db.ComDB.Close()
}

func openDbFile(add bool, dir, name string) (*os.File, error) {
	if !add {
		f, err := os.Create(path.Join(dir, name))
		if err != nil {
			return nil, err
		}
		return f, nil
	}

	f, err := os.OpenFile(path.Join(dir, name), os.O_RDWR, 0666)
	if err != nil {
		return nil, err
	}
	return f, nil
}
