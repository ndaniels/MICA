package cablastp

import (
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"os"
	"path"
	"strconv"
)

const (
	FileCoarseFasta = "coarse.fasta"
	FileCoarseLinks = "coarse.links"
	FileCoarseSeeds = "coarse.seeds"
	FileCompressed  = "compressed.cbp"
	FileIndex       = "index"
	FileParams      = "params"
)

type DBConf struct {
	MinMatchLen         int
	MatchKmerSize       int
	GappedWindowSize    int
	UngappedWindowSize  int
	ExtSeqIdThreshold   int
	MatchSeqIdThreshold int
	MatchExtend         int
	MapSeedSize         int
	ExtSeedSize         int
}

var DefaultDBConf = DBConf{
	MinMatchLen:         40,
	MatchKmerSize:       4,
	GappedWindowSize:    25,
	UngappedWindowSize:  10,
	ExtSeqIdThreshold:   50,
	MatchSeqIdThreshold: 60,
	MatchExtend:         30,
	MapSeedSize:         6,
	ExtSeedSize:         4,
}

func LoadDBConf(r io.Reader) (conf DBConf, err error) {
	defer func() {
		if perr := recover(); perr != nil {
			err = perr.(error)
		}
	}()
	conf = DefaultDBConf
	csvReader := csv.NewReader(r)
	csvReader.Comma = ':'
	csvReader.Comment = '#'
	csvReader.FieldsPerRecord = 2
	csvReader.TrailingComma = false
	csvReader.TrimLeadingSpace = true

	lines, err := csvReader.ReadAll()
	if err != nil {
		return conf, err
	}

	for _, line := range lines {
		atoi := func() int {
			var i64 int64
			var err error
			if i64, err = strconv.ParseInt(line[1], 10, 32); err != nil {
				panic(err)
			}
			return int(i64)
		}
		switch line[0] {
		case "MinMatchLen":
			conf.MinMatchLen = atoi()
		case "MatchKmerSize":
			conf.MatchKmerSize = atoi()
		case "GappedWindowSize":
			conf.GappedWindowSize = atoi()
		case "UngappedWindowSize":
			conf.UngappedWindowSize = atoi()
		case "ExtSeqIdThreshold":
			conf.ExtSeqIdThreshold = atoi()
		case "MatchSeqIdThreshold":
			conf.MatchSeqIdThreshold = atoi()
		case "MatchExtend":
			conf.MatchExtend = atoi()
		case "MapSeedSize":
			conf.MapSeedSize = atoi()
		case "ExtSeedSize":
			conf.ExtSeedSize = atoi()
		default:
			return conf, fmt.Errorf("Invalid DBConf flag: %s", line[0])
		}
	}

	return conf, nil
}

func (flagConf DBConf) FlagMerge(fileConf DBConf) (DBConf, error) {
	only := make(map[string]bool, 0)
	flag.Visit(func(f *flag.Flag) { only[f.Name] = true })

	if only["map-seed-size"] {
		return flagConf, fmt.Errorf("The map seed size cannot be changed for " +
			"an existing database.")
	}

	if !only["min-match-len"] {
		flagConf.MinMatchLen = fileConf.MinMatchLen
	}
	if !only["match-kmer-size"] {
		flagConf.MatchKmerSize = fileConf.MatchKmerSize
	}
	if !only["gapped-window-size"] {
		flagConf.GappedWindowSize = fileConf.GappedWindowSize
	}
	if !only["ungapped-window-size"] {
		flagConf.UngappedWindowSize = fileConf.UngappedWindowSize
	}
	if !only["ext-seq-id-threshold"] {
		flagConf.ExtSeqIdThreshold = fileConf.ExtSeqIdThreshold
	}
	if !only["match-seq-id-threshold"] {
		flagConf.MatchSeqIdThreshold = fileConf.MatchSeqIdThreshold
	}
	if !only["match-extend"] {
		flagConf.MatchExtend = fileConf.MatchExtend
	}
	if !only["ext-seed-size"] {
		flagConf.ExtSeedSize = fileConf.ExtSeedSize
	}
	return flagConf, nil
}

func (dbConf DBConf) Write(w io.Writer) error {
	csvWriter := csv.NewWriter(w)
	csvWriter.Comma = ':'
	csvWriter.UseCRLF = false

	s := func(i int) string {
		return fmt.Sprintf("%d", i)
	}
	records := [][]string{
		{"MinMatchLen", s(dbConf.MinMatchLen)},
		{"MatchKmerSize", s(dbConf.MatchKmerSize)},
		{"GappedWindowSize", s(dbConf.GappedWindowSize)},
		{"UngappedWindowSize", s(dbConf.UngappedWindowSize)},
		{"ExtSeqIdThreshold", s(dbConf.ExtSeqIdThreshold)},
		{"MatchSeqIdThreshold", s(dbConf.MatchSeqIdThreshold)},
		{"MatchExtend", s(dbConf.MatchExtend)},
		{"MapSeedSize", s(dbConf.MapSeedSize)},
		{"ExtSeedSize", s(dbConf.ExtSeedSize)},
	}
	if err := csvWriter.WriteAll(records); err != nil {
		return err
	}
	return nil
}

type DB struct {
	DBConf
	Name     string
	ComDB    *CompressedDB
	CoarseDB *CoarseDB

	coarseFasta, coarseSeeds, coarseLinks, compressed, index, params *os.File
}

func NewDB(conf DBConf, dir string) (*DB, error) {
	// Check to make sure the directory doesn't exist.
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

	err = os.Mkdir(dir, 0777)
	if err != nil {
		return nil, fmt.Errorf("Could not create directory '%s': %s.", dir, err)
	}

	db := &DB{
		DBConf: conf,
		Name:   path.Base(dir),
	}
	db.openDbFiles(dir, false)

	db.ComDB = NewCompressedDB(db.compressed, db.index)
	db.CoarseDB = NewCoarseDB(
		db.coarseFasta, db.coarseSeeds, db.coarseLinks, db.MapSeedSize)

	return db, nil
}

func NewAppendDB(conf DBConf, dir string) (*DB, error) {
	_, err := os.Open(dir)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s' for appending "+
			"because: %s.", dir, err)
	}

	db := &DB{
		DBConf: conf,
		Name:   path.Base(dir),
	}
	db.openDbFiles(dir, true)

	// Now try to load the configuration parameters from the 'params' file.
	// We always prefer options from 'params' except when it has been
	// overridden via the command line.
	paramConf, err := LoadDBConf(db.params)
	if err != nil {
		return nil, err
	}
	db.DBConf, err = db.DBConf.FlagMerge(paramConf)
	if err != nil {
		return nil, err
	}

	db.ComDB = NewAppendCompressedDB(db.compressed, db.index)
	db.CoarseDB = NewAppendCoarseDB(
		db.coarseFasta, db.coarseSeeds, db.coarseLinks, db.MapSeedSize)

	return db, nil
}

func LoadDB(dir string) (*DB, error) {
	_, err := os.Open(dir)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s' for reading "+
			"because: %s.", dir, err)
	}

	db := &DB{
		Name:        path.Base(dir),
		coarseSeeds: nil,
	}

	open := func(name string) (*os.File, error) {
		f, err := os.Open(path.Join(dir, name))
		if err != nil {
			return nil, err
		}
		return f, nil
	}

	db.coarseFasta, err = open(FileCoarseFasta)
	if err != nil {
		return nil, err
	}
	db.coarseLinks, err = open(FileCoarseLinks)
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
	db.params, err = open(FileParams)
	if err != nil {
		return nil, err
	}

	// Now try to load the configuration parameters from the 'params' file.
	// We always prefer options from 'params' except when it has been
	// overridden via the command line.
	db.DBConf, err = LoadDBConf(db.params)
	if err != nil {
		return nil, err
	}

	db.ComDB, err = LoadCompressedDB(db.compressed, db.index)
	if err != nil {
		return nil, err
	}

	db.CoarseDB, err = LoadCoarseDB(
		db.coarseFasta, db.coarseLinks, db.MapSeedSize)
	if err != nil {
		return nil, err
	}

	return db, nil
}

func (db *DB) Save() error {
	if err := db.CoarseDB.Save(); err != nil {
		return err
	}
	return nil
}

func (db *DB) ReadClose() {
	db.CoarseDB.ReadClose()
	db.ComDB.ReadClose()
}

func (db *DB) WriteClose() {
	db.CoarseDB.WriteClose()
	db.ComDB.WriteClose()
}

func (db *DB) openDbFiles(dir string, add bool) error {
	var err error

	// We need to open six files:
	// 1) Coarse FASTA database
	// 2) Coarse links
	// 3) Seeds
	// 4) Original links
	// 5) A line index
	// 6) A configuration file
	//
	// If 'add' is false, then each file is opened and truncated for writing.
	// If 'add' is true, then each file is opened in 'append' mode for writing.
	open := func(name string) (*os.File, error) {
		return openDbFile(add, dir, name)
	}

	db.coarseFasta, err = open(FileCoarseFasta)
	if err != nil {
		return err
	}
	db.coarseLinks, err = open(FileCoarseLinks)
	if err != nil {
		return err
	}
	db.coarseSeeds, err = open(FileCoarseSeeds)
	if err != nil {
		return err
	}
	db.compressed, err = open(FileCompressed)
	if err != nil {
		return err
	}
	db.index, err = open(FileIndex)
	if err != nil {
		return err
	}
	db.params, err = open(FileParams)
	if err != nil {
		return err
	}

	return nil
}

func openDbFile(add bool, dir, name string) (*os.File, error) {
	if !add {
		f, err := os.Create(path.Join(dir, name))
		if err != nil {
			return nil, err
		}
		return f, nil
	}

	f, err := os.OpenFile(path.Join(dir, name), os.O_RDWR|os.O_APPEND, 0666)
	if err != nil {
		return nil, err
	}
	return f, nil
}
