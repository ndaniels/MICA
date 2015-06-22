package main

import (
	"bytes"
	"fmt"
	"os"
	"runtime/pprof"

	"github.com/ndaniels/mica"
)

func s(i int) string {
	return fmt.Sprintf("%d", i)
}

func su(i uint64) string {
	return fmt.Sprintf("%d", i)
}

func writeFasta(oseqs []mica.OriginalSeq, buf *bytes.Buffer) error {
	for _, oseq := range oseqs {
		_, err := fmt.Fprintf(buf, "> %s\n%s\n",
			oseq.Name, string(oseq.Residues))
		if err != nil {
			return fmt.Errorf("Could not write to buffer: %s", err)
		}
	}
	return nil
}

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

func handleFatalError(msg string, err error) error {
	if err != nil {
		fatalf(msg+": %s\n", err)
	}
	return nil
}
