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

	newFastaIndex, err := os.Create("coarse.fasta.index.new")
		fatalf("Failed to create new index file: %s\n", err)
	}

	for i := 0; i < len(coarsedb.Seqs); i++ {
		buf.Reset()
		fmt.Fprintf(buf, ">%d\n%s\n", i, string(coarsedb.Seqs[i].Residues))

		err = binary.Write(newFastaIndex, binary.BigEndian, byteOff)
		if err != nil {
			fatalf("Failed to write to new index file: %s\n", err)
		}

		byteOff += int64(buf.Len())
	}

}

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}
