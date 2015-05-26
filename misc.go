package cablastp

import (
	"flag"
	"fmt"
	"os"
)

var (
	Verbose = false
)

func Vprint(s string) {
	if !Verbose {
		return
	}
	fmt.Fprint(os.Stderr, s)
}

func Vprintf(format string, v ...interface{}) {
	if !Verbose {
		return
	}
	fmt.Fprintf(os.Stderr, format, v...)
}

func Vprintln(s string) {
	if !Verbose {
		return
	}
	fmt.Fprintln(os.Stderr, s)
}

func PrintFlagDefaults() {
	flag.VisitAll(func(fg *flag.Flag) {
		fmt.Printf("--%s=\"%s\"\n\t%s\n", fg.Name, fg.DefValue, fg.Usage)
	})
}
