package cablastp

import (
	"fmt"
	"os"
	"path"
)

const (
	FileParams = "params"
)

type DB struct {
	DBConf
	Name     string
	Path     string
	ComDB    *CompressedDB
	CoarseDB *CoarseDB

	coarseFasta, coarseSeeds, coarseLinks, compressed, index, params *os.File
}

func NewWriteDB(appnd bool, conf DBConf, dir string) (*DB, error) {
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
		DBConf: conf,
		Name:   path.Base(dir),
		Path:   dir,
		params: nil,
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

	db.ComDB, err = NewWriteCompressedDB(appnd, db)
	if err != nil {
		return nil, err
	}
	db.CoarseDB, err = NewWriteCoarseDB(appnd, db)
	if err != nil {
		return nil, err
	}

	return db, nil
}

func (db *DB) openWriteFile(appnd bool, name string) (*os.File, error) {
	var f *os.File
	var err error

	if appnd {
		f, err = os.OpenFile(
			path.Join(db.Path, name), os.O_RDWR|os.O_APPEND, 0666)
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

func NewReadDB(dir string) (*DB, error) {
	_, err := os.Open(dir)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s' for reading "+
			"because: %s.", dir, err)
	}

	db := &DB{
		Name:        path.Base(dir),
		coarseSeeds: nil,
		params:      nil,
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
	db.ComDB, err = NewReadCompressedDB(db)
	if err != nil {
		return nil, err
	}
	db.CoarseDB, err = NewReadCoarseDB(db)
	if err != nil {
		return nil, err
	}

	return db, nil
}

func (db *DB) openReadFile(name string) (*os.File, error) {
	f, err := os.Open(path.Join(db.Path, name))
	if err != nil {
		return nil, err
	}
	return f, nil
}

func (db *DB) Save() error {
	if err := db.CoarseDB.Save(); err != nil {
		return err
	}
	if err := db.DBConf.Write(db.params); err != nil {
		return err
	}
	return nil
}

func (db *DB) ReadClose() {
	db.params.Close()
	db.CoarseDB.ReadClose()
	db.ComDB.ReadClose()
}

func (db *DB) WriteClose() {
	db.params.Close()
	db.CoarseDB.WriteClose()
	db.ComDB.WriteClose()
}
