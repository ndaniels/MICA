package main

import (
	"bytes"
	"encoding/xml"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path"
	"runtime"
	"runtime/pprof"

	// "code.google.com/p/biogo/io/seqio/fasta" 
	// "code.google.com/p/biogo/seq" 

	"github.com/BurntSushi/cablastp"
)

// blastp -db all-yeasts-cablastdb/blastdb -outfmt 5 < yeasts-query.fasta 

// A BLAST database is created on the reference sequence after compression.
// The search program will blast the input query sequences against this
// database (with a relaxed e-value), and expand the hits using links into an 
// in memory FASTA file. This FASTA file is passed to the stdin of a 
// `makeblastdb` command, which outputs the fine BLAST database. Finally, the 
// query sequences are blasted against this new database, and the hits are 
// returned unmodified.

const interval = 1000

var (
	// A default configuration.
	dbConf = cablastp.DefaultDBConf

	// Flags that affect the higher level operation of compression.
	// Flags that control algorithmic parameters are stored in `dbConf`.
	flagGoMaxProcs = runtime.NumCPU()
	flagQuiet      = false
	flagCpuProfile = ""
	flagMemProfile = ""
)

func init() {
	log.SetFlags(0)

	flag.StringVar(&dbConf.BlastMakeBlastDB, "makeblastdb",
		dbConf.BlastMakeBlastDB,
		"The location of the 'makeblastdb' executable.")
	flag.StringVar(&dbConf.BlastBlastp, "blastp",
		dbConf.BlastBlastp,
		"The location of the 'blastp' executable.")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagQuiet, "quiet", flagQuiet,
		"When set, the only outputs will be errors echoed to stderr.")
	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")

	flag.Usage = usage
	flag.Parse()

	runtime.GOMAXPROCS(flagGoMaxProcs)
}

func main() {
	buf := new(bytes.Buffer)

	if flag.NArg() < 2 {
		flag.Usage()
	}

	// If the quiet flag isn't set, enable verbose output.
	if !flagQuiet {
		cablastp.Verbose = true
	}

	queryFasta, err := os.Open(flag.Arg(1))
	if err != nil {
		fatalf("Could not open '%s': %s.\n", flag.Arg(1), err)
	}

	db, err := cablastp.NewReadDB(flag.Arg(0))
	if err != nil {
		fatalf("Could not open '%s' database: %s\n", flag.Arg(0), err)
	}
	cablastp.Vprintln("")

	cmd := exec.Command(
		db.BlastBlastp,
		"-db", path.Join(db.Path, cablastp.FileBlastCoarse),
		"-outfmt", "5")
	cmd.Stdout = buf
	cmd.Stdin = queryFasta
	if err := cablastp.Exec(cmd); err != nil {
		fatalf("%s\n", err)
	}

	results := blast{}
	if err := xml.NewDecoder(buf).Decode(&results); err != nil {
		fatalf("Could not parse BLAST search results: %s", err)
	}
	for _, hit := range results.Hits {
		for _, hsp := range hit.Hsps {
			oseqs, err := db.CoarseDB.Expand(db.ComDB,
				hit.Accession, hsp.HitFrom, hsp.HitTo)
			if err != nil {
				errorf("Could not decompress coarse sequence %d (%d, %d): %s\n",
					hit.Accession, hsp.HitFrom, hsp.HitTo, err)
				continue
			}

			fmt.Println("--------------------------------------")
			fmt.Printf("Accession: %d (%d, %d)\n",
				hit.Accession, hsp.HitFrom, hsp.HitTo)
			for _, oseq := range oseqs {
				fmt.Printf("%s\n\n", oseq)
			}
			fmt.Println("--------------------------------------")
		}
	}

	cleanup(db)
}

func cleanup(db *cablastp.DB) {
	if len(flagCpuProfile) > 0 {
		pprof.StopCPUProfile()
	}
	if len(flagMemProfile) > 0 {
		writeMemProfile(fmt.Sprintf("%s.last", flagMemProfile))
	}
	db.ReadClose()
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

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] "+
			"database-directory "+
			"query-fasta-file\n",
		path.Base(os.Args[0]))
	cablastp.PrintFlagDefaults()
	os.Exit(1)
}
