package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path"
	"runtime/pprof"
)

const (
	alphaSize = 20
)

var (
	flagInitDbLen          int
	flagMinMatchLen        int
	flagMatchKmerSize      int
	flagGappedWindowSize   int
	flagUngappedWindowSize int
	flagSeqIdThreshold     int
	flagSeedSize           int

	flagCpuProfile string
	flagMemProfile string
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&flagInitDbLen, "init-db-len", 1,
		"The initial size of the 'unique database' in terms of the number "+
			"of protein sequences added from the input database.")
	flag.IntVar(&flagMinMatchLen, "min-match-len", 25,
		"The minimum size of a match.")
	flag.IntVar(&flagMatchKmerSize, "match-kmer-size", 3,
		"The size of kmer fragments to match in ungapped extension.")
	flag.IntVar(&flagGappedWindowSize, "gapped-window-size", 10,
		"The size of the gapped match window.")
	flag.IntVar(&flagUngappedWindowSize, "ungapped-window-size", 10,
		"The size of the ungapped match window.")
	flag.IntVar(&flagSeqIdThreshold, "seq-id-threshold", 50,
		"The sequence identity threshold of a match. (An integer in the"+
			"inclusive range from 0 to 100.)")
	flag.IntVar(&flagSeedSize, "seed-size", 3,
		"The size of a seed.")

	flag.StringVar(&flagCpuProfile, "cpuprofile", "",
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", "",
		"When set, a memory profile will be written to the file specified.")
}

func main() {
	var err error

	if flag.NArg() < 1 {
		log.Println("At least one protein fasta file must be specified.")
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

	allseqs := make([][]*originalSeq, flag.NArg())
	for i, arg := range flag.Args() {
		allseqs[i], err = readSeqs(arg)
		if err != nil {
			log.Fatal(err)
		}
	}

	cdb := newCompressedDb([]*originalSeq{allseqs[0][0]})
	for _, origSeq := range allseqs[0][1:] {
		// "lastmatch" pointer to last residue of last match. (or 0)
		// "current" pointer to the starting residue of each kmer window. (or 0)

		// for every kmer, look up in seed table.
		// if no seed, increment "current" pointer by one.
		// else, start looking for a match.
		// ???
		// if no match is found, increment "current" pointer by one.
		// else, set "current" pointer to index of last residue + 1 in match.
		//		 and set "lastmatch" equal to "current" - 1.
		//		 and create link table entry ???

		// when do things get added to the compressed db?
		//
		// Boundary case: when the "current" pointer reaches end of current
		// sequence, any residues from it back to the "lastmatch" pointer are
		// added.
		//
		// After a match is found, all residues from "lastmatch" up to the
		// beginning of the match are added to the compressed db. (After
		// extending match in both directions.)
		//
		// "adding to the compressed db" includes:
		// 1. Adding a sequence to the set of sequences in the db.
		// 2. Adding all kmer windows in the sequence as seeds.
		cdb.add(origSeq)
	}

	fmt.Printf("%s\n", cdb)

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
		"Usage: %s [flags] fasta-file [fasta-file ...]",
		path.Base(os.Args[0]))
	flag.PrintDefaults()
	os.Exit(1)
}
