package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path"
	"runtime"
	"runtime/pprof"

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
	flagMakeBlastDB     = "makeblastdb"
	flagBlastp          = "blastp"
	flagDmnd            = "diamond"
	flagGoMaxProcs      = runtime.NumCPU()
	flagQuiet           = false
	flagCpuProfile      = ""
	flagMemProfile      = ""
	flagCoarseEval      = 5.0
	flagNoCleanup       = false
	flagDmndFine        = ""
	flagCoarseDmndMatch = 50
	flagFineDmndMatch   = 60
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

	flag.IntVar(&flagCoarseDmndMatch, "dmnd-coarse-match", flagCoarseDmndMatch,
		"The matching threshold for coarse search with diamond")
	flag.IntVar(&flagFineDmndMatch, "dmnd-fine-match", flagFineDmndMatch,
		"The matching threshold for fine search with diamond (assuming diamond fine search is enabled).")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")

	flag.IntVar(&flagCoarseDmndMatch, "dmnd-coarse-match", flagCoarseDmndMatch,
		"The matching threshold for coarse search with diamond")

	flag.IntVar(&flagFineDmndMatch, "dmnd-fine-match", flagFineDmndMatch,
		"The matching threshold for fine search with diamond (assuming diamond fine search is enabled).")

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

	db, err := neutronium.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}

	inputFastaQueryName := flag.Arg(1)
	aaQueryFile, err := os.Open(inputFastaQueryName)
	if err != nil {
		fatalf("Could not open '%s' query file: %s\n", inputFastaQueryName, err)
	}

	neutronium.Vprintln("\nProcessing Queries...")
	err = processQueries(db, aaQueryFile)
	if err != nil {
		fatalf("Error processing queries: %s\n", err)
	}

	cleanup(db)
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

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] database-directory query-fasta-file "+
			"[--blast-args BLASTP_ARGUMENTS]\n",
		path.Base(os.Args[0]))
	neutronium.PrintFlagDefaults()
	os.Exit(1)
}
