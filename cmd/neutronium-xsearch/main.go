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

	"github.com/TuftsBCB/seq"

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
	argDBConf = neutronium.DefaultDBConf.DeepCopy()
	// Flags that affect the operation of search.
	// Flags that control algorithmic parameters are stored in `queryDBConf`.
	flagMakeBlastDB    = "makeblastdb"
	flagBlastx         = "blastx"
	flagBlastn         = "blastn"
	flagDmnd           = "diamond"
	flagGoMaxProcs     = runtime.NumCPU()
	flagQuiet          = false
	flagCpuProfile     = ""
	flagMemProfile     = ""
	flagCoarseEval     = 5.0
	flagNoCleanup      = false
	flagCompressQuery  = false
	flagBatchQueries   = false
	flagIterativeQuery = false
)

// blastArgs are all the arguments after "--blast-args".
var blastArgs []string

func init() {
	log.SetFlags(0)

	// Regular cablastp-xsearch options

	flag.StringVar(&flagMakeBlastDB, "makeblastdb",
		flagMakeBlastDB,
		"The location of the 'makeblastdb' executable.")
	flag.StringVar(&flagBlastx, "blastx",
		flagBlastx,
		"The location of the 'blastx' executable.")
	flag.StringVar(&flagBlastn, "blastn",
		flagBlastn,
		"The location of the 'blastn' executable.")
	flag.StringVar(&flagDmnd, "diamond",
		flagDmnd,
		"The location of the 'diamond' executable.")
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
	flag.BoolVar(&flagIterativeQuery, "iterative-queries", flagIterativeQuery,
		"When set, will process queries one at a time instead of as a batch.")
	flag.BoolVar(&flagCompressQuery, "compress-query", flagCompressQuery,
		"When set, will process compress queries before search.")

	// compress options

	flag.IntVar(&argDBConf.MinMatchLen, "min-match-len",
		argDBConf.MinMatchLen,
		"The minimum size of a match.")
	flag.IntVar(&argDBConf.MatchKmerSize, "match-kmer-size",
		argDBConf.MatchKmerSize,
		"The size of kmer fragments to match in ungapped extension.")
	flag.IntVar(&argDBConf.ExtSeqIdThreshold, "ext-seq-id-threshold",
		argDBConf.ExtSeqIdThreshold,
		"The sequence identity threshold of [un]gapped extension. \n"+
			"\t(An integer in the inclusive range from 0 to 100.)")
	flag.IntVar(&argDBConf.MatchSeqIdThreshold, "match-seq-id-threshold",
		argDBConf.MatchSeqIdThreshold,
		"The sequence identity threshold of an entire match.")
	flag.IntVar(&argDBConf.MatchExtend, "match-extend",
		argDBConf.MatchExtend,
		"The maximum number of residues to blindly extend a \n"+
			"\tmatch without regard to sequence identity. This is \n"+
			"\tto avoid small sequences in the coarse database.")
	flag.IntVar(&argDBConf.MapSeedSize, "map-seed-size",
		argDBConf.MapSeedSize,
		"The size of a seed in the K-mer map. This size combined with\n"+
			"\t'ext-seed-size' forms the total seed size.")

	flag.IntVar(&argDBConf.LowComplexity, "low-complexity",
		argDBConf.LowComplexity,
		"The window size used to detect regions of low complexity.\n"+
			"\tLow complexity regions are repetitions of a single amino\n"+
			"\tacid residue. Low complexity regions are skipped when\n"+
			"\ttrying to extend a match.")
	flag.IntVar(&argDBConf.SeedLowComplexity, "seed-low-complexity",
		argDBConf.SeedLowComplexity,
		"The seed window size used to detect regions of low complexity.\n"+
			"\tLow complexity regions are repetitions of a single amino\n"+
			"\tacid residue. Low complexity regions matching this window\n"+
			"\tsize are not included in the seeds table.")
	// flag.Float64Var(&flagMaxSeedsGB, "max-seeds", flagMaxSeedsGB,
	// 	"When set, the in memory seeds table will be completely erased\n"+
	// 		"\twhen the memory used by seeds exceeds the specified number,\n"+
	// 		"\tin gigabytes.\n"+
	// 		"\tEach seed corresponds to 16 bytes of memory.\n"+
	// 		"\tSetting to zero disables this behavior.")

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

	db, err := neutronium.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}

	inputFastaQueryName := flag.Arg(1)
	nuclQueryFile, err := os.Open(inputFastaQueryName)
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", inputFastaQueryName, err)
	}

	neutronium.Vprintln("\nProcessing Queries...")
	err = processQueries(db, nuclQueryFile)
	if err != nil {
		fatalf("Error processing queries: %s\n", err)
	}

	cleanup(db)
}

func processQueries(db *neutronium.DB, nuclQueryFile *os.File) error {

	neutronium.Vprintln("\nBlasting with diamond query on coarse database...")
	dmndOutFile, err := dmndBlastXCoarse(db, nuclQueryFile)
	if err != nil {
		fatalf("Error blasting with diamond on coarse database: %s\n", err)
	}

	neutronium.Vprintln("Decompressing diamond hits...")
	dmndOutArr, err := ioutil.ReadAll(dmndOutFile)
	if err != nil {
		return fmt.Errorf("Could not read diamond output: %s", err)
	}
	if len(dmndOutArr) == 0 {
		return fmt.Errorf("No coarse hits. %s", "Aborting.")
	}
	neutronium.Vprintln("Expanding diamond hits...")
	dmndOut := bytes.NewBuffer(dmndOutArr)
	expandedSequences, err := expandDmndHits(db, dmndOut)
	if err != nil {
		fatalf("%s\n", err)
	}

	// Write the contents of the expanded sequences to a fasta file.
	// It is then indexed using makeblastdb.
	searchBuf := new(bytes.Buffer)
	if err := writeFasta(expandedSequences, searchBuf); err != nil {
		fatalf("Could not create FASTA input from coarse hits: %s\n", err)
	}

	// Create the fine blast db in a temporary directory
	neutronium.Vprintln("Building fine BLAST database...")
	tmpFineDB, err := makeFineBlastDB(db, searchBuf)
	handleFatalError("Could not create fine database to search on", err)

	// retrieve the cluster members for the original representative query seq

	// pass them to blastx on the expanded (fine) db

	// Finally, run the query against the fine fasta database and pass on the
	// stdout and stderr...
	bs, err := ioutil.ReadAll(nuclQueryFile)
	if err != nil {
		return fmt.Errorf("Could not read input fasta query: %s", err)
	}
	nuclQueryReader := bytes.NewReader(bs)

	err = blastFine(db, tmpFineDB, nuclQueryReader)
	handleFatalError("Error blasting fine database", err)

	// Delete the temporary fine database.
	if !flagNoCleanup {
		err := os.RemoveAll(tmpFineDB)
		handleFatalError("Could not delete fine BLAST database", err)
	}
	return nil
}

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

	cmd := exec.Command(flagBlastx, flags...)
	cmd.Stdin = stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
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

func expandCoarseSequence(db *neutronium.DB, seqId int, coarseSequence *seq.Sequence) ([]neutronium.OriginalSeq, error) {
	originalSeqs, err := db.CoarseDB.Expand(db.ComDB, seqId, 0, coarseSequence.Len())
	if err != nil {
		return nil, err
	}

	return originalSeqs, nil
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

func dmndBlastXCoarse(db *neutronium.DB, queries *os.File) (*os.File, error) {
	// diamond blastp -d nr -q reads.fna -a matches -t <temporary directory>

	dmndOutFile, err := ioutil.TempFile(".", "dmnd-out-")
	if err != nil {
		return nil, fmt.Errorf("Could not build temporary file for diamond output: %s", err)
	}

	cmd := exec.Command(
		flagDmnd,
		"blastx",
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

func getInputFasta(inputFilename string) (*bytes.Reader, error) {
	queryFasta, err := os.Open(inputFilename)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s': %s.", flag.Arg(1), err)
	}
	bs, err := ioutil.ReadAll(queryFasta)
	if err != nil {
		return nil, fmt.Errorf("Could not read input fasta query: %s", err)
	}
	return bytes.NewReader(bs), nil
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

func handleFatalError(msg string, err error) error {
	if err != nil {
		fatalf(msg+": %s\n", err)
	}
	return nil
}
