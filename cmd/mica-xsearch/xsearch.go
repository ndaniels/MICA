package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"io"
	"os"
	"compress/gzip"
	"strings"
	"github.com/ndaniels/mica"
)

func processQueries(db *mica.DB, inputFastaQueryName string) error {


	mica.Vprintln("\nBlasting with diamond query on coarse database...")

	dmndOutDaaFilename, err := dmndBlastXCoarse(db, inputFastaQueryName)
	if err != nil {
		return fmt.Errorf("Error blasting with diamond on coarse database: %s\n", err)
	}


	dmndOutFile, err := convertDmndToBlastTabular(dmndOutDaaFilename)
	if err != nil {
		return fmt.Errorf("Error converting diamond output to blast tabular: %s\n", err)
	}

	mica.Vprintln("Decompressing diamond hits...")
	dmndOutArr, err := ioutil.ReadAll(dmndOutFile)

	if !flagNoCleanup {
		err := os.RemoveAll(dmndOutFile.Name())
		handleFatalError("Could not delete diamond output from coarse search", err)
		err = os.RemoveAll(dmndOutDaaFilename)
		handleFatalError("Could not delete diamond output from coarse search", err)
	}

	if err != nil {
		return fmt.Errorf("Could not read diamond output: %s", err)
	}
	if len(dmndOutArr) == 0 {
		return fmt.Errorf("No coarse hits. %s", "Aborting.")
	}
	mica.Vprintln("Expanding diamond hits...")
	dmndOut := bytes.NewBuffer(dmndOutArr)
	expandedSequences, err := expandDmndHits(db, dmndOut)
	if err != nil {
		return fmt.Errorf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	searchBuf := new(bytes.Buffer)
	if err := writeFasta(expandedSequences, searchBuf); err != nil {
		fatalf("Could not create FASTA input from coarse hits: %s\n", err)
	}

	if flagDmndFine != "" {

		mica.Vprintln("Building fine DIAMOND database...")
		tmpFineDB, err := makeFineDmndDB(searchBuf)
		handleFatalError("Could not create fine diamond database to search on", err)

		err = dmndBlastXFine(inputFastaQueryName, flagDmndFine, tmpFineDB)
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

		var nuclQueryFile io.Reader
		nuclQueryFile, err = os.Open(inputFastaQueryName)
		if err != nil {
			fatalf("Could not open '%s' query file: %s\n", inputFastaQueryName, err)
		}

		if strings.HasSuffix(inputFastaQueryName, ".gz") {
			nuclQueryFile, err = gzip.NewReader(nuclQueryFile)
			if err != nil {
				return err
			}
		}

		// Finally, run the query against the fine fasta database and pass on the
		// stdout and stderr...
		bs, err := ioutil.ReadAll(nuclQueryFile)
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
