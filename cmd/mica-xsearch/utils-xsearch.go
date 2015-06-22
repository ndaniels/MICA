package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"runtime/pprof"

	"github.com/TuftsBCB/seq"

	"github.com/ndaniels/mica"
)

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}

func errorf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func writeMemProfile(name string) {
	f, err := os.Create(name)
	if err != nil {
		fatalf("%s\n", err)
	}
	pprof.WriteHeapProfile(f)
	f.Close()
}

func expandCoarseSequence(db *mica.DB, seqId int, coarseSequence *seq.Sequence) ([]mica.OriginalSeq, error) {
	originalSeqs, err := db.CoarseDB.Expand(db.ComDB, seqId, 0, coarseSequence.Len())
	if err != nil {
		return nil, err
	}

	return originalSeqs, nil
}

func getInputFasta(inputFilename string) (*bytes.Reader, error) {
	queryFasta, err := os.Open(inputFilename)
	if err != nil {
		return nil, fmt.Errorf("Could not open '%s': %s.", flag.Arg(1), err)
	}
	bs, err := ioutil.ReadAll(queryFasta)
	if err != nil {
		return nil, fmt.Errorf("Could not read input fasta query: %s", err)
	}
	return bytes.NewReader(bs), nil
}

func s(i int) string {
	return fmt.Sprintf("%d", i)
}

func su(i uint64) string {
	return fmt.Sprintf("%d", i)
}

func writeFasta(oseqs []mica.OriginalSeq, buf *bytes.Buffer) error {
	for _, oseq := range oseqs {
		_, err := fmt.Fprintf(buf, ">%s\n%s\n",
			oseq.Name, string(oseq.Residues))
		if err != nil {
			return fmt.Errorf("Could not write to buffer: %s", err)
		}
	}
	return nil
}

func handleFatalError(msg string, err error) error {
	if err != nil {
		fatalf(msg+": %s\n", err)
	}
	return nil
}
