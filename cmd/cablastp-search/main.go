package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path"
	"runtime"
	"runtime/pprof"

	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq"

	"github.com/BurntSushi/cablastp"
)

// blastp -db all-yeasts-cablastdb/blastdb -outfmt 5 < yeasts-query.fasta 

// A BLAST database is created on the reference sequence after compression.
// The search program will blast the input query sequences against this
// database (with a relaxed e-value), and expand the hits using links into an 
// in memory FASTA file. This FASTA file is passed to the stdin of a 
// `makeblastdb` command, which outputs the fine BLAST database. Finally, the 
// query sequences are blasted against this new database, and the hits are 
// returned unmodified.

const interval = 1000

var (
	// Flags that affect the higher level operation of compression.
	// Flags that control algorithmic parameters are stored in `dbConf`.
	flagGoMaxProcs = runtime.NumCPU()
	flagQuiet      = false
	flagCpuProfile = ""
	flagMemProfile = ""
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagQuiet, "quiet", flagQuiet,
		"When set, the only outputs will be errors echoed to stderr.")
	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")

	flag.Usage = usage
	flag.Parse()

	runtime.GOMAXPROCS(flagGoMaxProcs)
}

func main() {
	if flag.NArg() < 2 {
		flag.Usage()
	}

	// If the quiet flag isn't set, enable verbose output.
	if !flagQuiet {
		cablastp.Verbose = true
	}

	// Open the fasta file specified for writing.
	// This will be the source of the fine database.
	outFasta, err := os.Create(flag.Arg(1))
	if err != nil {
		fatalf("Could not write to '%s': %s\n", flag.Arg(1), err)
	}
	fastaWriter := fasta.NewWriter(outFasta, 60)

	// Create a new database for writing. If we're appending, we load
	// the coarse database into memory, and setup the database for writing.
	db, err := cablastp.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}
	cablastp.Vprintln("")

	// Start the CPU profile after all of the data has been read.
	if len(flagCpuProfile) > 0 {
		f, err := os.Create(flagCpuProfile)
		if err != nil {
			fatalf("%s\n", err)
		}
		pprof.StartCPUProfile(f)
	}

	numSeqs := db.ComDB.NumSequences()
	for orgSeqId := 0; orgSeqId < numSeqs; orgSeqId++ {
		oseq, err := db.ComDB.ReadNextSeq(db.CoarseDB, orgSeqId)
		if err != nil {
			fatalf("%s\n", err)
		}
		fastaWriter.Write(
			seq.New(oseq.Name, append(oseq.Residues, '*'), nil))
	}

	cleanup(db)
	if err = fastaWriter.Close(); err != nil {
		fatalf("%s\n", err)
	}
}

func cleanup(db *cablastp.DB) {
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
		"\nUsage: %s [flags] "+
			"database-directory "+
			"out-fasta-file\n",
		path.Base(os.Args[0]))
	cablastp.PrintFlagDefaults()
	os.Exit(1)
}
