package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path"
	"runtime/pprof"

	"github.com/BurntSushi/cablastp"
)

var (
	flagMinMatchLen        = 25
	flagMatchKmerSize      = 3
	flagGappedWindowSize   = 25
	flagUngappedWindowSize = 10
	flagSeqIdThreshold     = 50
	flagSeedSize           = 3

	flagCpuProfile = ""
	flagMemProfile = ""
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&flagMinMatchLen, "min-match-len", flagMinMatchLen,
		"The minimum size of a match.")
	flag.IntVar(&flagMatchKmerSize, "match-kmer-size", flagMatchKmerSize,
		"The size of kmer fragments to match in ungapped extension.")
	flag.IntVar(&flagGappedWindowSize, "gapped-window-size",
		flagGappedWindowSize, "The size of the gapped match window.")
	flag.IntVar(&flagUngappedWindowSize, "ungapped-window-size",
		flagUngappedWindowSize, "The size of the ungapped match window.")
	flag.IntVar(&flagSeqIdThreshold, "seq-id-threshold", flagSeqIdThreshold,
		"The sequence identity threshold of a match. (An integer in the "+
			"inclusive range from 0 to 100.)")
	flag.IntVar(&flagSeedSize, "seed-size", flagSeedSize,
		"The size of a seed.")

	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")
}

func main() {
	var err error

	if flag.NArg() < 3 {
		flag.Usage()
	}

	if len(flagCpuProfile) > 0 {
		f, err := os.Create(flagCpuProfile)
		if err != nil {
			fatalf("%s\n", err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// For each FASTA file, convert each sequence in each FASTA file to
	// an OriginalSeq. This preps them for processing and conversion into
	// CompressedSeq.
	orgSeqId := 0
	refdb := newReferenceDB()
	comdb := cablastp.NewCompressedDB()
	for _, arg := range flag.Args()[3:] {
		err = cablastp.ReadOriginalSeqs(arg,
			func(orgSeq *cablastp.OriginalSeq) {
				comSeq := compress(refdb, orgSeqId, orgSeq)
				comdb.Add(comSeq)
				orgSeqId++
				if orgSeqId % 100 == 0 {
					fmt.Printf("%d complete\n", orgSeqId)
				}
			})
		if err != nil {
			log.Fatal(err)
		}
	}

	// Initialize the reference and compressed databases. For each original
	// sequence, convert it to a compressed sequence and add it to the
	// compressed database. (The process of compressing an original sequence
	// will add to the reference database if applicable.)
	// orgSeqId := 0 
	// for _, fastaSeqs := range allseqs { 
		// for _, orgSeq := range fastaSeqs { 
			// comSeq := compress(refdb, orgSeqId, orgSeq) 
			// comdb.Add(comSeq) 
			// orgSeqId++ 
			// if orgSeqId % 100 == 0 { 
				// fmt.Printf("%d complete\n", orgSeqId) 
			// } 
		// } 
	// } 

	if err := refdb.savePlain(flag.Arg(0), flag.Arg(1)); err != nil {
		fatalf("Could not save coarse database: %s\n", err)
	}
	if err := comdb.SavePlain(flag.Arg(2)); err != nil {
		fatalf("Could not save compressed database: %s\n", err)
	}

	if len(flagMemProfile) > 0 {
		f, err := os.Create(flagMemProfile)
		if err != nil {
			fatalf("%s\n", err)
		}
		pprof.WriteHeapProfile(f)
		f.Close()
	}
}

func errorf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func fatalf(format string, v ...interface{}) {
	errorf(format, v...)
	os.Exit(1)
}

func init() {
	flag.Usage = usage
	flag.Parse()
}

func usage() {
	fmt.Fprintf(os.Stderr,
		"Usage: %s [flags] "+
			"output-coarse-fasta-name "+
			"output-compressed-links "+
			"fasta-file [fasta-file ...]\n",
		path.Base(os.Args[0]))
	flag.PrintDefaults()
	os.Exit(1)
}
