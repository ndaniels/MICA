package main

import (
	"bytes"
	"encoding/xml"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path"

	"github.com/ndaniels/neutronium"
)

func expandBlastHits(
	db *neutronium.DB, blastOut *bytes.Buffer) ([]neutronium.OriginalSeq, error) {

	results := blast{}
	if err := xml.NewDecoder(blastOut).Decode(&results); err != nil {
		return nil, fmt.Errorf("Could not parse BLAST search results: %s", err)
	}

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]neutronium.OriginalSeq, 0, 100)
	for _, hit := range results.Hits {

		for _, hsp := range hit.Hsps {
			someOseqs, err := db.CoarseDB.Expand(db.ComDB,
				hit.Accession, hsp.HitFrom, hsp.HitTo)
			if err != nil {
				errorf("Could not decompress coarse sequence %d (%d, %d): %s\n",
					hit.Accession, hsp.HitFrom, hsp.HitTo, err)
				continue
			}

			// Make sure this hit is below the coarse e-value threshold.
			if hsp.Evalue > flagCoarseEval {
				continue
			}

			for _, oseq := range someOseqs {
				if used[oseq.Id] {
					continue
				}
				used[oseq.Id] = true
				oseqs = append(oseqs, oseq)
			}
		}
	}
	// if len(oseqs) == 0 {
	// 	return nil, fmt.Errorf("No hits from coarse search\n")
	// }
	return oseqs, nil
}

func blastCoarse(
	db *neutronium.DB, stdin *bytes.Reader, stdout *bytes.Buffer) error {

	cmd := exec.Command(
		flagBlastn,
		"-db", path.Join(db.Path, neutronium.FileBlastCoarse),
		"-num_threads", s(flagGoMaxProcs),
		"-outfmt", "5", "-dbsize", su(db.BlastDBSize))
	cmd.Stdin = stdin
	cmd.Stdout = stdout
	return neutronium.Exec(cmd)
}

func makeFineBlastDB(db *neutronium.DB, stdin *bytes.Buffer) (string, error) {
	tmpDir, err := ioutil.TempDir("", "cablastp-fine-search-db")
	if err != nil {
		return "", fmt.Errorf("Could not create temporary directory: %s\n", err)
	}

	cmd := exec.Command(
		flagMakeBlastDB, "-dbtype", "prot",
		"-title", neutronium.FileBlastFine,
		"-in", "-",
		"-out", path.Join(tmpDir, neutronium.FileBlastFine))
	cmd.Stdin = stdin

	neutronium.Vprintf("Created temporary fine BLAST database in %s\n", tmpDir)

	return tmpDir, neutronium.Exec(cmd)
}

func blastFine(
	db *neutronium.DB, blastFineDir string, stdin *bytes.Reader) error {

	// We pass our own "-db" flag to blastp, but the rest come from user
	// defined flags.
	flags := []string{
		"-db", path.Join(blastFineDir, neutronium.FileBlastFine),
		"-dbsize", su(db.BlastDBSize),
		"-num_threads", s(flagGoMaxProcs),
	}
	flags = append(flags, blastArgs...)

	cmd := exec.Command(flagBlastx, flags...)
	cmd.Stdin = stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return neutronium.Exec(cmd)
}
