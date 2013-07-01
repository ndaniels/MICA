package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path"
	"runtime"
	"runtime/pprof"

	"github.com/TuftsBCB/io/fasta"

	"github.com/BurntSushi/cablastp"
)

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
	outFasta, err := os.Create(flag.Arg(1))
	if err != nil {
		fatalf("Could not write to '%s': %s\n", flag.Arg(1), err)
	}
	fastaWriter := fasta.NewWriter(outFasta)
	fastaWriter.Asterisk = true

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
		oseq, err := db.ComDB.ReadSeq(db.CoarseDB, orgSeqId)
		if err != nil {
			fatalf("Error reading seq id '%d': %s\n", orgSeqId, err)
		}
		if err := fastaWriter.Write(oseq.FastaSeq()); err != nil {
			cablastp.Vprintf("Error writing seq '%s': %s\n", oseq.Name, err)
		}
	}

	cleanup(db)
	if err = fastaWriter.Flush(); err != nil {
		fatalf("%s\n", err)
	}
	if err = outFasta.Close(); err != nil {
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
