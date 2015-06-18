package main

import (
	"bytes"
	"encoding/binary"
	"flag"
	"fmt"
	"github.com/ndaniels/neutronium"
	"os"
)

func init() {

	flag.Parse()
}

func main() {
	db, err := neutronium.NewReadDB(flag.Arg(0))

	byteOff := int64(0)
	buf := new(bytes.Buffer)
	coarsedb := db.CoarseDB

	newFastaFile, err := os.Create("coarse.fasta.new")
	if err != nil {
		fatalf("Failed to create new sequence file: %s\n", err)
	}

	newFastaIndex, err := os.Create("coarse.fasta.index.new")
	if err != nil {
		fatalf("Failed to create new index file: %s\n", err)
	}

	for i, seq := range coarsedb.Seqs {
		neutronium.Vprintln(seq.String())
		buf.Reset()
		fmt.Fprintf(buf, ">%d\n%s\n", i, seq.String())
		if _, err = newFastaFile.Write(buf.Bytes()); err != nil {
			return
		}

		err = binary.Write(newFastaIndex, binary.BigEndian, byteOff)
		if err != nil {
			fatalf("Failed to write to new index file: %s\n", err)
		}

		byteOff += int64(buf.Len())
	}

	newFastaFile.Close()
	newFastaIndex.Close()

}

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}
