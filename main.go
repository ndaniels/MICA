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
	flagInitDbLen      int
	flagMinMatchLen    int
	flagSeqIdThreshold float64
	flagSeedSize int
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&flagInitDbLen, "init-db-len", 10000,
		"The initial size of the 'unique database'.")
	flag.IntVar(&flagMinMatchLen, "min-match-len", 300,
		"The minimum size of a match.")
	flag.Float64Var(&flagSeqIdThreshold, "seq-id-threshold", 0.8,
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
