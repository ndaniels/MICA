package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/ndaniels/neutronium"
)

func expandDmndHits(db *neutronium.DB, dmndOut *bytes.Buffer) ([]neutronium.OriginalSeq, error) {

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]neutronium.OriginalSeq, 0, 100)

	dmndScanner := bufio.NewScanner(dmndOut)
	for dmndScanner.Scan() {
		line := dmndScanner.Text()
		if err := dmndScanner.Err(); err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}

		// Example line:
		// 0        1          2             3          4              5             6           7         8             9           10    11
		// queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore
		// YAL001C  897745     96.12         1160       45             0             1           1160      1             1160        0e+00 2179.8

		splitLine := strings.Fields(line)

		if len(splitLine) < 12 {
			return nil, fmt.Errorf("Line in diamond output is too short: %s", line)
		}

		coarseID, err := strconv.Atoi(splitLine[1])
		if err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		hitFrom, err := strconv.Atoi(splitLine[8])
		if err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		hitTo, err := strconv.Atoi(splitLine[9])
		if err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		eval, err := strconv.ParseFloat(splitLine[10], 64)
		if err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}

		someOseqs, err := db.CoarseDB.Expand(db.ComDB, coarseID, hitFrom, hitTo)
		if err != nil {
			return nil, fmt.Errorf("Could not decompress coarse sequence %d (%d, %d): %s\n", coarseID, hitFrom, hitTo, err)
		}

		// Make sure this hit is below the coarse e-value threshold.
		if eval > flagCoarseEval {
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
	return oseqs, nil
}

func dmndCoarse(db *neutronium.DB, queries *os.File) (*os.File, error) {
	// diamond blastp -d nr -q reads.fna -a matches -t <temporary directory>

	dmndOutFile, err := ioutil.TempFile(".", "dmnd-out-")
	if err != nil {
		return nil, fmt.Errorf("Could not build temporary file for diamond output: %s", err)
	}

	cmd := exec.Command(
		flagDmnd,
		"blastp",
		"-d", path.Join(db.Path, neutronium.FileDmndCoarse),
		"-q", queries.Name(),
		"--threads", s(flagGoMaxProcs),
		"-o", dmndOutFile.Name(),
		"--compress", "0",
		"--top", "90")
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout

	err = neutronium.Exec(cmd)
	if err != nil {
		return nil, fmt.Errorf("Error using diamond to blast coarse db: %s", err)
	}

	return dmndOutFile, nil
}
