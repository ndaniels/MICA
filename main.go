package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
)

const (
	alphaSize = 20
)

var (
	flagInitDbLen         int
	flagMinMatchLen       int
	flagMatchKmerSize int
	flagUngappedWindowSize int
	flagSeqIdThreshold    float64
	flagSeedSize          int
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&flagInitDbLen, "init-db-len", 10000,
		"The initial size of the 'unique database'.")
	flag.IntVar(&flagMinMatchLen, "min-match-len", 100,
		"The minimum size of a match.")
	flag.IntVar(&flagMatchKmerSize, "match-kmer-size", 2,
		"The size of kmer fragments to match in ungapped extension.")
	flag.IntVar(&flagUngappedWindowSize, "ungapped-window-size", 4,
		"The size of the ungapped match window.")
	flag.Float64Var(&flagSeqIdThreshold, "seq-id-threshold", 0.5,
		"The sequence identity threshold of a match")
	flag.IntVar(&flagSeedSize, "seed-size", 3,
		"The size of a seed.")
}

func main() {
	var err error

	flag.Usage = usage
	flag.Parse()
	if flag.NArg() < 1 {
		log.Println("At least one protein fasta file must be specified.")
		flag.Usage()
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
	fmt.Println("")
	fmt.Printf("%s\n", cdb.seeds)
}

func usage() {
	basename := os.Args[0]
	if lastSlash := strings.LastIndex(basename, "/"); lastSlash > -1 {
		basename = basename[lastSlash+1:]
	}
	log.Printf("Usage: %s [flags] fasta-file [fasta-file ...]", basename)
	flag.PrintDefaults()
	os.Exit(1)
}
