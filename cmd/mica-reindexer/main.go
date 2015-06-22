package main

import (
	"bytes"
	"encoding/binary"
	"flag"
	"fmt"
	"github.com/ndaniels/mica"
	"os"
)

func init() {

	flag.Parse()
}

func main() {
	db, err := mica.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Failed to open database: %s\n", err)
	}

	byteOff := int64(0)
	buf := new(bytes.Buffer)
	newFastaFile, err := os.Create("coarse.fasta.new")
	if err != nil {
		fatalf("Failed to create new sequence file: %s\n", err)
	}

	newFastaIndex, err := os.Create("coarse.fasta.index.new")
	if err != nil {
		fatalf("Failed to create new index file: %s\n", err)
	}

	err = db.CoarseDB.LoadSeqs()
	if err != nil {
		fatalf("Failed to load coarse db sequences into memory: %s\n", err)
	}

	for i, seq := range db.CoarseDB.Seqs {
		buf.Reset()
		fmt.Fprintf(buf, ">%d\n%s\n", i, string(seq.Residues))
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
