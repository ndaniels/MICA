package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"

	"github.com/ndaniels/neutronium"
)

func processQueries(db *neutronium.DB, inputFastaQueryFile *os.File) error {
	neutronium.Vprintln("\nBlasting with diamond query on coarse database...")
	dmndOutFile, err := dmndCoarse(db, inputFastaQueryFile)
	if err != nil {
		return fmt.Errorf("Error blasting with diamond on coarse database: %s\n", err)
	}

	neutronium.Vprintln("Decompressing diamond hits...")
	dmndOutArr, err := ioutil.ReadAll(dmndOutFile)
	if err != nil {
		return fmt.Errorf("Could not read diamond output: %s\n", err)
	}
	if len(dmndOutArr) == 0 {
		return fmt.Errorf("No coarse hits: %s\n", "Aborting.")
	}
	neutronium.Vprintln("Expanding diamond hits...")
	dmndOut := bytes.NewBuffer(dmndOutArr)
	expandedSequences, err := expandDmndHits(db, dmndOut)
	if err != nil {
		return fmt.Errorf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	buf := new(bytes.Buffer)
	if err := writeFasta(expandedSequences, buf); err != nil {
		return fmt.Errorf("Could not create FASTA input from coarse hits: %s\n", err)
	}

	// Create the fine blast db in a temporary directory
	neutronium.Vprintln("Building fine BLAST database...")
	tmpDir, err := makeFineBlastDB(db, buf)
	if err != nil {
		return fmt.Errorf("Could not create fine database to search on: %s\n", err)
	}

	// Finally, run the query against the fine fasta database and pass on the
	// stdout and stderr...
	inputFastaQuery, err := getInputFasta()
	if err != nil {
		return fmt.Errorf("Could not read input fasta query: %s\n", err)
	}
	neutronium.Vprintln("Blasting query on fine database...")
	if _, err := inputFastaQuery.Seek(0, os.SEEK_SET); err != nil {
		return fmt.Errorf("Could not seek to start of query fasta input: %s\n", err)
	}
	if err := blastFine(db, tmpDir, inputFastaQuery); err != nil {
		return fmt.Errorf("Error blasting fine database: %s\n", err)
	}

	// Delete the temporary fine database.
	if !flagNoCleanup {
		if err := os.RemoveAll(dmndOutFile.Name()); err != nil {
			return fmt.Errorf("Could not delete output from diamond: %s\n", err)
		}
		if err := os.RemoveAll(tmpDir); err != nil {
			return fmt.Errorf("Could not delete fine BLAST database: %s\n", err)
		}
	}

	return nil
}
