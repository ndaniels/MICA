package main

import (
	"bufio"
	"bytes"
	"encoding/xml"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"

	"github.com/ndaniels/neutronium"
)

// A BLAST database is created on the reference sequence after compression.
// The search program will blast the input query sequences against this
// database (with a relaxed e-value), and expand the hits using links into an
// in memory FASTA file. This FASTA file is passed to the stdin of a
// `makeblastdb` command, which outputs the fine BLAST database. Finally, the
// query sequences are blasted against this new database, and the hits are
// returned unmodified.

var (
	// A default configuration.
	dbConf = neutronium.DefaultDBConf

	// Flags that affect the operation of search.
	// Flags that control algorithmic parameters are stored in `dbConf`.
	flagMakeBlastDB = "makeblastdb"
	flagBlastp      = "blastp"
	flagDmnd        = "diamond"
	flagGoMaxProcs  = runtime.NumCPU()
	flagQuiet       = false
	flagCpuProfile  = ""
	flagMemProfile  = ""
	flagCoarseEval  = 5.0
	flagNoCleanup   = false
)

// blastArgs are all the arguments after "--blast-args".
var blastArgs []string

func init() {
	log.SetFlags(0)

	flag.StringVar(&flagMakeBlastDB, "makeblastdb",
		flagMakeBlastDB,
		"The location of the 'makeblastdb' executable.")
	flag.StringVar(&flagBlastp, "blastp",
		flagBlastp,
		"The location of the 'blastp' executable.")
	flag.StringVar(&flagDmnd, "dmnd-blastp",
		flagDmnd,
		"The location of the 'diamond-blastp' executable.")
	flag.Float64Var(&flagCoarseEval, "coarse-eval", flagCoarseEval,
		"The e-value threshold for the coarse search. This will NOT\n"+
			"\tbe used on the fine search. The fine search e-value threshold\n"+
			"\tcan be set in the 'blast-args' argument.")
	flag.BoolVar(&flagNoCleanup, "no-cleanup", flagNoCleanup,
		"When set, the temporary fine BLAST database that is created\n"+
			"\twill NOT be deleted.")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagQuiet, "quiet", flagQuiet,
		"When set, the only outputs will be errors echoed to stderr.")
	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")

	// find '--blast-args' and chop off the remainder before letting the flag
	// package have its way.
	for i, arg := range os.Args {
		if arg == "--blast-args" {
			blastArgs = os.Args[i+1:]
			os.Args = os.Args[:i]
		}
	}

	flag.Usage = usage
	flag.Parse()

	runtime.GOMAXPROCS(flagGoMaxProcs)
}

func main() {

	if flag.NArg() != 2 {
		flag.Usage()
	}

	// If the quiet flag isn't set, enable verbose output.
	if !flagQuiet {
		neutronium.Verbose = true
	}

	inputFastaQueryFile, err := getInputFastaFile()
	if err != nil {
		fatalf("Could not open input fasta query: %s\n", err)
	}

	db, err := neutronium.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}

	neutronium.Vprintln("\nBlasting with diamond query on coarse database...")
	dmndOutFile, err := dmndCoarse(db, inputFastaQueryFile)
	if err != nil {
		fatalf("Error blasting with diamond on coarse database: %s\n", err)
	}

	neutronium.Vprintln("Decompressing diamond hits...")
	dmndOutArr, err := ioutil.ReadAll(dmndOutFile)
	if err != nil {
		fatalf("Could not read diamond output: %s", err)
	}
	dmndOut := bytes.NewBuffer(dmndOutArr)
	expandedSequences, err := expandDmndHits(db, dmndOut)
	if err != nil {
		fatalf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	buf := new(bytes.Buffer)
	if err := writeFasta(expandedSequences, buf); err != nil {
		fatalf("Could not create FASTA input from coarse hits: %s\n", err)
	}

	// Create the fine blast db in a temporary directory
	neutronium.Vprintln("Building fine BLAST database...")
	tmpDir, err := makeFineBlastDB(db, buf)
	if err != nil {
		fatalf("Could not create fine database to search on: %s\n", err)
	}

	// Finally, run the query against the fine fasta database and pass on the
	// stdout and stderr...
	inputFastaQuery, err := getInputFasta()
	if err != nil {
		fatalf("Could not read input fasta query: %s\n", err)
	}
	neutronium.Vprintln("Blasting query on fine database...")
	if _, err := inputFastaQuery.Seek(0, os.SEEK_SET); err != nil {
		fatalf("Could not seek to start of query fasta input: %s\n", err)
	}
	if err := blastFine(db, tmpDir, inputFastaQuery); err != nil {
		fatalf("Error blasting fine database: %s\n", err)
	}

	// Delete the temporary fine database.
	if !flagNoCleanup {
		if err := os.RemoveAll(tmpDir); err != nil {
			fatalf("Could not delete fine BLAST database: %s\n", err)
		}
	}

	cleanup(db)
}

func s(i int) string {
	return fmt.Sprintf("%d", i)
}

func su(i uint64) string {
	return fmt.Sprintf("%d", i)
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

	cmd := exec.Command(flagBlastp, flags...)
	cmd.Stdin = stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return neutronium.Exec(cmd)
}

func makeFineBlastDB(db *neutronium.DB, stdin *bytes.Buffer) (string, error) {
	tmpDir, err := ioutil.TempDir("", "neutronium-fine-search-db")
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

func writeFasta(oseqs []neutronium.OriginalSeq, buf *bytes.Buffer) error {
	for _, oseq := range oseqs {
		_, err := fmt.Fprintf(buf, "> %s\n%s\n",
			oseq.Name, string(oseq.Residues))
		if err != nil {
			return fmt.Errorf("Could not write to buffer: %s", err)
		}
	}
	return nil
}

func expandDmndHits(db *neutronium.DB, dmndOut *bytes.Buffer) ([]neutronium.OriginalSeq, error) {
	scanner := bufio.NewScanner(os.Stdin)
	for scanner.Scan() {
		fmt.Println(scanner.Text()) // Println will add back the final '\n'
	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, "reading standard input:", err)
	}

	used := make(map[int]bool, 100) // prevent original sequence duplicates
	oseqs := make([]neutronium.OriginalSeq, 0, 100)

	dmndScanner := bufio.NewScanner(dmndOut)
	for dmndScanner.Scan() {
		line := dmndScanner.Text()
		if err := dmndScanner.Err(); err != nil {
			return nil, fmt.Errorf("Error reading from diamond output: %s", err)
		}

		coarseID := -1
		hitFrom := -1
		hitTo := -1
		eval := -1.0

		for i, word := range strings.Split(line, "\t") {
			if i == 1 {
				_coarseID, err := strconv.Atoi(word)
				if err != nil {
					return nil, fmt.Errorf("Error reading from diamond output: %s", err)
				}
				coarseID = _coarseID
			} else if i == 8 {
				_hitFrom, err := strconv.Atoi(word)
				if err != nil {
					return nil, fmt.Errorf("Error reading from diamond output: %s", err)
				}
				hitFrom = _hitFrom
			} else if i == 9 {
				_hitTo, err := strconv.Atoi(word)
				if err != nil {
					return nil, fmt.Errorf("Error reading from diamond output: %s", err)
				}
				hitTo = _hitTo
			} else if i == 10 {
				_eval, err := strconv.ParseFloat(word, 64)
				if err != nil {
					return nil, fmt.Errorf("Error reading from diamond output: %s", err)
				}
				eval = _eval
			}
		}

		someOseqs, err := db.CoarseDB.Expand(db.ComDB, coarseID, hitFrom, hitTo)
		if err != nil {
			errorf("Could not decompress coarse sequence %d (%d, %d): %s\n", coarseID, hitFrom, hitTo, err)
			continue
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

func expandBlastHits(db *neutronium.DB, blastOut *bytes.Buffer) ([]neutronium.OriginalSeq, error) {

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
	return oseqs, nil
}

func blastCoarse(
	db *neutronium.DB, stdin *bytes.Reader, stdout *bytes.Buffer) error {

	cmd := exec.Command(
		flagBlastp,
		"-db", path.Join(db.Path, neutronium.FileBlastCoarse),
		"-num_threads", s(flagGoMaxProcs),
		"-outfmt", "5", "-dbsize", su(db.BlastDBSize))
	cmd.Stdin = stdin
	cmd.Stdout = stdout
	return neutronium.Exec(cmd)
}

func dmndCoarse(db *neutronium.DB, queries *os.File) (*os.File, error) {
	// diamond blastx -d nr -q reads.fna -a matches -t <temporary directory>

	dmndOutFile, err := ioutil.TempFile(".", "dmnd-out-")
	if err != nil {
		return nil, fmt.Errorf("Could not build temporary file for diamond output: %s", err)
	}

	// dmndCmd := fmt.Sprintf("%s blastp -d %s -q %s --threads %d -o %s ",
	// 						flagDmnd,
	// 						path.Join(db.Path, neutronium.FileDmndCoarse),
	// 						queries.Name(),
	// 						s(flagGoMaxProcs),
	// 						)

	cmd := exec.Command(
		flagDmnd,
		"blastp",
		"-d", path.Join(db.Path, neutronium.FileDmndCoarse),
		"-q", queries.Name(),
		"--threads", s(flagGoMaxProcs),
		"-o", dmndOutFile.Name())
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout

	err = neutronium.Exec(cmd)
	if err != nil {
		return nil, fmt.Errorf("Error using diamond to blast coarse db: %s", err)
	}

	return dmndOutFile, nil
}

func getInputFasta() (*bytes.Reader, error) {
	queryFasta, err := os.Open(flag.Arg(1))
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s': %s.", flag.Arg(1), err)
	}
	bs, err := ioutil.ReadAll(queryFasta)
	if err != nil {
		return nil, fmt.Errorf("Could not read input fasta query: %s", err)
	}
	return bytes.NewReader(bs), nil
}

func getInputFastaFile() (*os.File, error) {
	queryFasta, err := os.Open(flag.Arg(1))
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s': %s.", flag.Arg(1), err)
	}
	return queryFasta, nil
}

func cleanup(db *neutronium.DB) {
	if len(flagCpuProfile) > 0 {
		pprof.StopCPUProfile()
	}
	if len(flagMemProfile) > 0 {
		writeMemProfile(fmt.Sprintf("%s.last", flagMemProfile))
	}
	db.ReadClose()
}

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}

func errorf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func writeMemProfile(name string) {
	f, err := os.Create(name)
	if err != nil {
		fatalf("%s\n", err)
	}
	pprof.WriteHeapProfile(f)
	f.Close()
}

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] database-directory query-fasta-file "+
			"[--blast-args BLASTP_ARGUMENTS]\n",
		path.Base(os.Args[0]))
	neutronium.PrintFlagDefaults()
	os.Exit(1)
}
