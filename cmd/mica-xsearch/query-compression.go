package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"

	"github.com/ndaniels/mica"
)

func processCompressedQueries(db *mica.DB, nuclQueryFileLoc string) error {

	queryDbLoc, err := compressQueries(nuclQueryFileLoc)
	if err != nil {
		return fmt.Errorf("Error compressing queries: %s\n",err)
	}
	
	queryDb, err := mica.NewReadDB(queryDbLoc)
	if err != nil {
		return fmt.Errorf("Error opening newly created compressed query db: %s\n", err)
	}


	mica.Vprintln("\nBlasting with diamond query on coarse database...")
	dmndOutDaaFilename, err := dmndBlastXCoarse(db, queryDb.CoarseDB.FileFasta.Name())
	if err != nil {
		return fmt.Errorf("Error blasting with diamond on coarse database from compressed queries: %s\n", err)
	}
	


	dmndOutFile, err := convertDmndToBlastTabular(dmndOutDaaFilename)
	if err != nil {
		return fmt.Errorf("Error converting diamond output to blast tabular: %s\n")
	}

	mica.Vprintln("Decompressing diamond hits...")
	dmndOutArr, err := ioutil.ReadAll(dmndOutFile)


	// nuclQueryFile, err := os.Open(nuclQueryFileLoc)
	// if err != nil {
	// 	fatalf("Could not open '%s' query file for fine search: %s\n", nuclQueryFileLoc, err)
	// }

	if !flagNoCleanup {
		err := os.RemoveAll(dmndOutFile.Name())
		handleFatalError("Could not delete diamond output from coarse search", err)
		err = os.RemoveAll(dmndOutDaaFilename)
		handleFatalError("Could not delete diamond output from coarse search", err)
		err = os.RemoveAll(queryDbLoc)
		handleFatalError("Could not delete diamond output from coarse search", err)
	}



	if err != nil {
		return fmt.Errorf("Could not read diamond output: %s", err)
	}
	if len(dmndOutArr) == 0 {
		return fmt.Errorf("No coarse hits. %s", "Aborting.")
	}
	mica.Vprintln("Expanding diamond hits (queries and targets)...")
	dmndOut := bytes.NewBuffer(dmndOutArr)
	expandedSequences, expandedQueries, err := expandDmndHitsAndQuery(db, queryDb, dmndOut)
	if err != nil {
		return fmt.Errorf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	searchBuf := new(bytes.Buffer)
	if err := writeFasta(expandedSequences, searchBuf); err != nil {
		fatalf("Could not create FASTA input from coarse hits: %s\n", err)
	}


	qSearchBuf := new(bytes.Buffer)
	if err := writeFasta(expandedQueries, qSearchBuf); err != nil {
		fatalf("Could not create FASTA input from coarse QUERY hits: %s\n", err)
	}
	fineQueryFile, err := ioutil.TempFile(flagTempFileDir, "fine-query-sequences")
	if err != nil {
		fatalf("Could not create temporary QUERY sequence file: %s\n", err)
	}
	err = ioutil.WriteFile(fineQueryFile.Name(), qSearchBuf.Bytes(), 0666)
	if err != nil {
		fatalf("Could not write to temporary QUERY sequence file: %s\n", err)
	}

	if flagDmndFine != "" {

		mica.Vprintln("Building fine DIAMOND database...")
		tmpFineDB, err := makeFineDmndDB(searchBuf)
		handleFatalError("Could not create fine diamond database to search on", err)

		err = dmndBlastXFine(fineQueryFile.Name(), flagDmndFine, tmpFineDB)
		handleFatalError("Error diamond-blasting (x-search) fine database", err)

		// Delete the temporary fine database.
		if !flagNoCleanup {
			err := os.RemoveAll(tmpFineDB)
			err = os.RemoveAll(tmpFineDB + ".dmnd")
			handleFatalError("Could not delete fine DIAMOND database", err)
		}

	} else {

		// Create the fine blast db in a temporary directory
		mica.Vprintln("Building fine BLAST database...")
		tmpFineDB, err := makeFineBlastDB(db, searchBuf)
		handleFatalError("Could not create fine blast database to search on", err)

		// retrieve the cluster members for the original representative query seq

		// pass them to blastx on the expanded (fine) db

		// Finally, run the query against the fine fasta database and pass on the
		// stdout and stderr...
		bs, err := ioutil.ReadAll(fineQueryFile)
		if err != nil {
			return fmt.Errorf("Could not read input fasta query: %s", err)
		}
		nuclQueryReader := bytes.NewReader(bs)

		err = blastFine(db, tmpFineDB, nuclQueryReader)
		handleFatalError("Error blasting fine database (x-search):", err)

		// Delete the temporary fine database.
		if !flagNoCleanup {
			err := os.RemoveAll(tmpFineDB)
			handleFatalError("Could not delete fine BLAST database", err)
		}
	}

	return nil
}


func compressQueries(queryFileName string) (string, error) {
	// dbDirLoc, err := ioutil.TempDir(flagTempFileDir, "temporary-compressed-query-db")
	dbDirLoc := flagTempFileDir 
	dbDirLoc = dbDirLoc + "/temporary-compressed-query-db"

	if flagQueryDBConf != "" {
		qdbParams, err := os.Open(flagQueryDBConf)
		if err != nil {
			return "", fmt.Errorf("Failed to load query db conf: %s", err)
		}
		queryDBConf, err = mica.LoadDBConf(qdbParams)
		if err != nil {
			return "", fmt.Errorf("Failed to load query db conf: %s", err)
		}
	}

	db, err := mica.NewWriteDB(false, queryDBConf, dbDirLoc)
	handleFatalError("Failed to open new db", err)
	mica.Vprintln("Starting query compress workers...")
	pool := mica.StartCompressReducedWorkers(db)
	seqId := db.ComDB.NumSequences()
	mainQuit := make(chan struct{}, 0)

	seqChan, err := mica.ReadOriginalSeqs(queryFileName, []byte{})
	handleFatalError("Could not read query sequences", err)
	mica.Vprintln("Reading sequences into query database...")
	for readSeq := range seqChan {
		// Do a non-blocking receive to see if main needs to quit.
		select {
		case <-mainQuit:
			<-mainQuit // wait for cleanup to finish before exiting main.
			return "", nil
		default:
		}

		handleFatalError("Failed to read sequence", readSeq.Err)

		queryDBConf.BlastDBSize += uint64(readSeq.Seq.Len())
		redReadSeq := &mica.ReducedSeq{
			&mica.Sequence{
				Name:     readSeq.Seq.Name,
				Residues: readSeq.Seq.Residues,
				Offset:   readSeq.Seq.Offset,
				Id:       readSeq.Seq.Id,
			},
		}
		seqId = pool.CompressReduced(seqId, redReadSeq)
	}
	mica.Vprintln("Cleaning up query database...")
	mica.CleanupDB(db, &pool)
	mica.Vprintln("")

	return dbDirLoc, nil
}