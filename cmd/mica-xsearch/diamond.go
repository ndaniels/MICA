package main

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/ndaniels/mica"
	"io/ioutil"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"
)

func dmndBlastXFine(queryFilename string, outFilename, fineFilename string) error {
     	cmd := exec.Command(
		flagDmnd,
		"blastx",
		"--sensitive",
		"-d", fineFilename,
		"-q", queryFilename,
		"--threads", s(flagGoMaxProcs),
		"-a", outFilename,
		"--compress", "0",
		"--top", s(flagFineDmndMatch))
		
     	if flagFast {     
	   cmd = exec.Command(
		flagDmnd,
		"blastx",
		"-d", fineFilename,
		"-q", queryFilename,
		"--threads", s(flagGoMaxProcs),
		"-a", outFilename,
		"--compress", "0",
		"--top", s(flagFineDmndMatch))
	}
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout

	err := mica.Exec(cmd)
	if err != nil {
		return fmt.Errorf("Error using diamond to blast coarse db: %s\n", err)
	}
	if !flagDmndOutput {
		daaFilename := outFilename + ".daa"
		tabularFile, err := convertDmndToBlastTabular(daaFilename)
		if err != nil {
			return fmt.Errorf("Error converting diamond output: %s\n", err)
		}
		os.Rename(tabularFile.Name(), outFilename)

	}

	return nil
}

func dmndBlastXCoarse(db *mica.DB, queryFilename string) (string, error) {
	// diamond blastp -d nr -q reads.fna -a matches -t <temporary directory>


	dmndOutFile, err := ioutil.TempFile(flagTempFileDir, "dmnd-coarse-out-")
	if err != nil {
		return "", fmt.Errorf("Could not open file to catch coarse search output: %s", err)
	}
	dmndOutFilename := dmndOutFile.Name()



	cmd := exec.Command(
		flagDmnd,
		"blastx",
		"--sensitive",
		"-d", path.Join(db.Path, mica.FileDmndCoarse),
		"-q", queryFilename,
		"--threads", s(flagGoMaxProcs),
		"-a", dmndOutFilename,
		"--compress", "0",
		"--top", s(flagCoarseDmndMatch),
		"--tmpdir", flagTempFileDir)

	if flagFast{
	   cmd = exec.Command(
		flagDmnd,
		"blastx",
		"-d", path.Join(db.Path, mica.FileDmndCoarse),
		"-q", queryFilename,
		"--threads", s(flagGoMaxProcs),
		"-a", dmndOutFilename,
		"--compress", "0",
		"--top", s(flagCoarseDmndMatch),
		"--tmpdir", flagTempFileDir)
	}
	
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout

	errp := mica.Exec(cmd)
	if errp != nil {
		return "", fmt.Errorf("Error using diamond to blast coarse db: %s", errp)
	}

	dmndOutFilename = dmndOutFilename + ".daa"
	if _, err := os.Stat(dmndOutFilename); err != nil {
	   return "", fmt.Errorf("Diamond did not produce daa output file: %s (%s)", errp, err)
	}
	
	return dmndOutFilename, nil
}

func convertDmndToBlastTabular(daaName string) (*os.File, error) {
	dmndOutFile, err := ioutil.TempFile(flagTempFileDir, "dmnd-out-tab-")
	if err != nil {
		return nil, fmt.Errorf("Could not build temporary file for diamond output: %s", err)
	}


	cmd := exec.Command(
		flagDmnd,
		"view",
		"-o", dmndOutFile.Name(),
		"-a", daaName)

	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout

	err = mica.Exec(cmd)
	if err != nil {
		return nil, fmt.Errorf("Error converting daa file to blast tabular: %s", err)
	}
	// err = os.Remove(daa.Name())
	// if err != nil {
	// 	return nil, fmt.Errorf("Error destroying .daa file: %s", err)
	// }
	return dmndOutFile, nil
}

func expandDmndHits(db *mica.DB, dmndOut *bytes.Buffer) ([]mica.OriginalSeq, error) {

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]mica.OriginalSeq, 0, 100)
	i := 0
	j := 0 
	dmndScanner := bufio.NewScanner(dmndOut)
	if flagPreload{
	   db.CoarseDB.Preload()
	   }
	for dmndScanner.Scan() {
	    	if i % 1000 == 0 {
		 mica.Vprintf("\rScanning %d Writing %d...",i,j)
	    	   
		   }
		i = i + 1
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
		       if j % 1000 == 0 {
	    	       	  mica.Vprintf("\rScanning %d Writing %d...",i,j)
		   	}
			j = j + 1
			if used[oseq.Id] {
				continue
			}
			used[oseq.Id] = true
			oseqs = append(oseqs, oseq)
		}
	}
	if flagPreload{
	   db.CoarseDB.Unload()
	   }
	return oseqs, nil
}

func expandDmndHitsAndQuery(db *mica.DB, qdb *mica.DB, dmndOut *bytes.Buffer) ([]mica.OriginalSeq, []mica.OriginalSeq, error) {

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]mica.OriginalSeq, 0, 100)
	qUsed := make(map[int]bool, 100)
	qSeqs := make([]mica.OriginalSeq, 0, 100)

	dmndScanner := bufio.NewScanner(dmndOut)
	
	for dmndScanner.Scan() {
		line := dmndScanner.Text()
		if err := dmndScanner.Err(); err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}

		// Example line:
		// 0        1          2             3          4              5             6           7         8             9           10    11
		// queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore
		// YAL001C  897745     96.12         1160       45             0             1           1160      1             1160        0e+00 2179.8

		splitLine := strings.Fields(line)

		if len(splitLine) < 12 {
			return nil, nil, fmt.Errorf("Line in diamond output is too short: %s", line)
		}
		coarseQID, err := strconv.Atoi(splitLine[0])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		coarseID, err := strconv.Atoi(splitLine[1])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		qHitFrom, err := strconv.Atoi(splitLine[6])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		qHitTo, err := strconv.Atoi(splitLine[7])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		hitFrom, err := strconv.Atoi(splitLine[8])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		hitTo, err := strconv.Atoi(splitLine[9])
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}
		eval, err := strconv.ParseFloat(splitLine[10], 64)
		if err != nil {
			return nil, nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}

		someOseqs, err := db.CoarseDB.Expand(db.ComDB, coarseID, hitFrom, hitTo)
		if err != nil {
			return nil, nil, fmt.Errorf("Could not decompress coarse sequence %d (%d, %d): %s\n", coarseID, hitFrom, hitTo, err)
		}

		someQseqs, err := qdb.CoarseDB.Expand(qdb.ComDB, coarseQID, qHitFrom, qHitTo)
		if err != nil {
			return nil, nil, fmt.Errorf("Could not decompress coarse  QUERY sequence %d (%d, %d): %s\n", coarseQID, qHitFrom, qHitTo, err)
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

		for _, qseq := range someQseqs {
			if qUsed[qseq.Id] {
				continue
			}
			qUsed[qseq.Id] = true
			qSeqs = append(qSeqs, qseq)
		}
	}
	return oseqs, qSeqs, nil
}

func makeFineDmndDB(seqBuf *bytes.Buffer) (string, error) {
	tmpSeqFile, err := ioutil.TempFile(flagTempFileDir, "fine-sequences-")
	if err != nil {
		return "", fmt.Errorf("Could not create temporary sequence file: %s\n", err)
	}
	err = ioutil.WriteFile(tmpSeqFile.Name(), seqBuf.Bytes(), 0666)
	if err != nil {
		return "", fmt.Errorf("Could not write to temporary sequence file: %s\n", err)
	}
	tmpDmndFile, err := ioutil.TempFile(flagTempFileDir, "fine-dmnd-db-")
	if err != nil {
		return "", fmt.Errorf("Could not create temporary diamond file: %s\n", err)
	}
	cmd := exec.Command(
		flagDmnd,
		"makedb",
		"--in", tmpSeqFile.Name(),
		"-d", tmpDmndFile.Name())

	err = mica.Exec(cmd)
	if err != nil {
		return "", fmt.Errorf("Could not create fine diamond database: %s\n", err)
	}

	err = os.RemoveAll(tmpSeqFile.Name())
	if err != nil {
		return "", fmt.Errorf("Could not remove temporary sequence file: %s\n", err)
	}

	return tmpDmndFile.Name(), nil

}

