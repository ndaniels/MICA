package cablastp

import (
	"fmt"
	"os"
	"os/exec"
	"path"
	"strings"
)

const (
	FileParams      = "params"
	FileBlastCoarse = "blastdb-coarse"
	FileBlastFine   = "blastdb-fine"
)

// A DB represents a cablastp database, which has three main components:
// a coarse database, a compressed database and a configuration file.
//
// A DB can be opened either for writing/appending (compression) or for
// reading (decompression).
type DB struct {
	// An embedded configuration.
	*DBConf

	// The path to the directory on disk.
	Path string

	// The name of the database. This corresponds to the basename of the path.
	Name string

	// The compressed database component.
	ComDB *CompressedDB

	// The coarse database component.
	CoarseDB *CoarseDB

	// Whether the database is being appended to. This can only occur if the
	// database was *initially* created without the read-only flag set.
	appending bool

	// File pointers.
	coarseFasta, coarseSeeds, coarseLinks, compressed, index, params *os.File
}

// NewWriteDB creates a new cablastp database, and prepares it for writing (or
// opens an existing database and prepares it for appending if 'appnd' is set).
//
// An error is returned if there is a problem accessing any of the files in
// the database.
//
// It is an error to open a database for writing that already exists if
// 'appnd' is not set.
//
// 'conf' should be a database configuration, typically defined (initially) from
// command line parameters. Note that if 'appnd' is set, then the configuration
// will be read from disk---only options explicitly set via the command line
// will be overwritten.
func NewWriteDB(appnd bool, conf *DBConf, dir string) (*DB, error) {
	Vprintf("Opening database in %s...\n", dir)

	if strings.HasSuffix(dir, ".tar") || strings.HasSuffix(dir, ".gz") {
		return nil, fmt.Errorf("The CaBLASTP database you've provided does " +
			"not appear to be a directory. Please make sure you've extracted " +
			"the downloaded database with `tar zxf cablastp-xxx.tar.gz` " +
			"before using it with CaBLASTP.")
	}

	_, err := os.Open(dir)
	if appnd {
		if err != nil {
			return nil, fmt.Errorf("Could not open '%s' for appending "+
				"because: %s.", dir, err)
		}
	} else {
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

		err = os.Mkdir(dir, 0777)
		if err != nil {
			return nil,
				fmt.Errorf("Could not create directory '%s': %s.", dir, err)
		}
	}

	db := &DB{
		DBConf:    conf,
		Name:      path.Base(dir),
		Path:      dir,
		params:    nil,
		appending: appnd,
	}

	// Do a sanity check and make sure we can access the `makeblastdb`
	// executable. Otherwise we might do a lot of work for nothing...
	if err = execExists(db.BlastMakeBlastDB); err != nil {
		return nil, fmt.Errorf(
			"Could not find 'makeblastdb' executable: %s", err)
	}

	// Now try to load the configuration parameters from the 'params' file.
	// We always prefer options from 'params' except when it has been
	// overridden via the command line.
	// We only need to do this when appending, otherwise a 'params' file
	// does not exist yet.
	db.params, err = db.openWriteFile(appnd, FileParams)
	if err != nil {
		return nil, err
	}
	if appnd {
		// If we're appending, we need some way of merging the configuration
		// on disk and the configuration supplied by the user. The hack here is
		// to inspect the command line arguments, and override the existing
		// configuration only with options explicitly set on the command line.
		paramConf, err := LoadDBConf(db.params)
		if err != nil {
			return nil, err
		}
		db.DBConf, err = db.DBConf.FlagMerge(paramConf)
		if err != nil {
			return nil, err
		}

		// If it's a read only database, we can't append!
		if db.ReadOnly {
			return nil, fmt.Errorf("Appending to a read-only database is " +
				"not possible.")
		}
	}

	db.ComDB, err = newWriteCompressedDB(appnd, db)
	if err != nil {
		return nil, err
	}
	db.CoarseDB, err = newWriteCoarseDB(appnd, db)
	if err != nil {
		return nil, err
	}

	Vprintf("Done opening database in %s.\n", dir)
	return db, nil
}

func (db *DB) filePath(name string) string {
	return path.Join(db.Path, name)
}

func (db *DB) openWriteFile(appnd bool, name string) (*os.File, error) {
	var f *os.File
	var err error

	if appnd {
		f, err = os.OpenFile(
			path.Join(db.Path, name), os.O_RDWR, 0666)
		if err != nil {
			return nil, err
		}
	} else {
		f, err = os.Create(path.Join(db.Path, name))
		if err != nil {
			return nil, err
		}
	}
	return f, nil
}

// NewReadDB opens a cablastp database for reading. An error is returned if
// there is a problem accessing any of the files on disk.
//
// Also, if the 'makeblastdb' or 'blastp' executales are not found, then an
// error is returned.
func NewReadDB(dir string) (*DB, error) {
	Vprintf("Opening database in %s...\n", dir)

	if strings.HasSuffix(dir, ".tar") || strings.HasSuffix(dir, ".gz") {
		return nil, fmt.Errorf("The CaBLASTP database you've provided does " +
			"not appear to be a directory. Please make sure you've extracted " +
			"the downloaded database with `tar zxf cablastp-xxx.tar.gz` " +
			"before using it with CaBLASTP.")
	}

	_, err := os.Open(dir)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s' for reading "+
			"because: %s.", dir, err)
	}

	db := &DB{
		Name:        path.Base(dir),
		Path:        dir,
		coarseSeeds: nil,
		params:      nil,
		appending:   false,
	}

	db.params, err = db.openReadFile(FileParams)
	if err != nil {
		return nil, err
	}

	// Now try to load the configuration parameters from the 'params' file.
	db.DBConf, err = LoadDBConf(db.params)
	if err != nil {
		return nil, err
	}

	// Do a sanity check and make sure we can access the `makeblastdb`
	// and `blastp` executables. Otherwise we might do a lot of work for
	// nothing...
	if err = execExists(db.BlastMakeBlastDB); err != nil {
		return nil, fmt.Errorf(
			"Could not find 'makeblastdb' executable: %s", err)
	}

	db.ComDB, err = newReadCompressedDB(db)
	if err != nil {
		return nil, err
	}
	db.CoarseDB, err = newReadCoarseDB(db)
	if err != nil {
		return nil, err
	}

	Vprintf("Done opening database in %s.\n", dir)
	return db, nil
}

func (db *DB) openReadFile(name string) (*os.File, error) {
	f, err := os.Open(path.Join(db.Path, name))
	if err != nil {
		return nil, err
	}
	return f, nil
}

// Save will write the contents of the database to disk. This should be called
// after compression is complete.
//
// After the database is saved, a blastp database is created from the coarse
// database.
//
// N.B. The compressed database is written as each sequence is processed, so
// this call will only save the coarse database. This may take a *very* long
// time if the database is not read only (since the seeds table has to be
// written).
func (db *DB) Save() error {
	var err error

	// Only write the params file when the database is first created.
	if !db.appending {
		// Make sure the params file is truncated so that we overwrite any
		// previous configuration.
		if err = db.params.Truncate(0); err != nil {
			return err
		}
		if _, err = db.params.Seek(0, os.SEEK_SET); err != nil {
			return err
		}
		if err = db.DBConf.Write(db.params); err != nil {
			return err
		}
	}

	// Write the coarse database to disk.
	// We don't need to explicitly save the compressed database, since its
	// data is written as it is generated (including the index).
	if err = db.CoarseDB.save(); err != nil {
		return err
	}

	// Now we need to construct a blastp database from the coarse fasta file.
	// e.g., `makeblastdb -dbtype prot -in coarse.fasta`
	cmd := exec.Command(
		db.BlastMakeBlastDB, "-dbtype", "prot",
		"-in", path.Join(db.Path, FileCoarseFasta),
		"-out", path.Join(db.Path, FileBlastCoarse))

	Vprintf("Creating %s...\n", FileBlastCoarse)
	if err = Exec(cmd); err != nil {
		return err
	}
	Vprintf("Done creating %s.\n", FileBlastCoarse)
	return nil
}

// ReadClose closes all appropriate files after reading from a database.
func (db *DB) ReadClose() {
	db.params.Close()
	db.CoarseDB.readClose()
	db.ComDB.readClose()
}

// WriteClose closes all appropriate files after writing to a database.
func (db *DB) WriteClose() {
	db.params.Close()
	db.CoarseDB.writeClose()
	db.ComDB.writeClose()
}

// execExists tests whether a binary exists in one's PATH.
func execExists(name string) error {
	_, err := exec.LookPath(name)
	if err != nil {
		return err
	}
	return nil
}

func fileExists(name string) error {
	_, err := os.Stat(name)
	if err != nil && os.IsNotExist(err) {
		return err
	}
	return nil
}
