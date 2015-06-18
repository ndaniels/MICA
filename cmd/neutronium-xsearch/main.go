package main

import (
	"flag"
	"fmt"
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
	argDBConf = neutronium.DefaultDBConf.DeepCopy()
	// Flags that affect the operation of search.
	// Flags that control algorithmic parameters are stored in `queryDBConf`.
	flagMakeBlastDB     = "makeblastdb"
	flagBlastx          = "blastx"
	flagBlastn          = "blastn"
	flagDmnd            = "diamond"
	flagGoMaxProcs      = runtime.NumCPU()
	flagQuiet           = false
	flagCpuProfile      = ""
	flagMemProfile      = ""
	flagCoarseEval      = 5.0
	flagNoCleanup       = false
	flagCompressQuery   = false
	flagBatchQueries    = false
	flagIterativeQuery  = false
	flagDmndFine        = ""
	flagCoarseDmndMatch = 90
	flagFineDmndMatch   = 90
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
	flag.StringVar(&flagDmndFine, "dmnd-fine",
		"",
		"When set, will use diamond for fine search writing the results to the specified file")

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
