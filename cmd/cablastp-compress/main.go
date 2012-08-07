package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"os/signal"
	"path"
	"runtime"
	"runtime/pprof"
	"time"

	"github.com/BurntSushi/cablastp"
)

const interval = 1000

var (
	timer           time.Time
	ignoredResidues = []byte{'J', 'O', 'U'}
	dbConf          = cablastp.DefaultDBConf

	flagGoMaxProcs  = runtime.NumCPU()
	flagAppend      = false
	flagOverwrite   = false
	flagVerbose     = false
	flagCpuProfile  = ""
	flagMemProfile  = ""
	flagMemStats    = ""
	flagMemInterval = false
)

func init() {
	log.SetFlags(0)

	flag.IntVar(&dbConf.MinMatchLen, "min-match-len",
		dbConf.MinMatchLen,
		"The minimum size of a match.")
	flag.IntVar(&dbConf.MatchKmerSize, "match-kmer-size",
		dbConf.MatchKmerSize,
		"The size of kmer fragments to match in ungapped extension.")
	flag.IntVar(&dbConf.GappedWindowSize, "gapped-window-size",
		dbConf.GappedWindowSize,
		"The size of the gapped match window.")
	flag.IntVar(&dbConf.UngappedWindowSize, "ungapped-window-size",
		dbConf.UngappedWindowSize,
		"The size of the ungapped match window.")
	flag.IntVar(&dbConf.ExtSeqIdThreshold, "ext-seq-id-threshold",
		dbConf.ExtSeqIdThreshold,
		"The sequence identity threshold of [un]gapped extension. \n"+
			"\t(An integer in the inclusive range from 0 to 100.)")
	flag.IntVar(&dbConf.MatchSeqIdThreshold, "match-seq-id-threshold",
		dbConf.MatchSeqIdThreshold,
		"The sequence identity threshold of an entire match.")
	flag.IntVar(&dbConf.MatchExtend, "match-extend",
		dbConf.MatchExtend,
		"The maximum number of residues to blindly extend a \n"+
			"\tmatch without regard to sequence identity. This is \n"+
			"\tto avoid small sequences in the coarse database.")
	flag.IntVar(&dbConf.MapSeedSize, "map-seed-size",
		dbConf.MapSeedSize,
		"The size of a seed in the K-mer map. This size combined with\n"+
			"\t'ext-seed-size' forms the total seed size.")
	flag.IntVar(&dbConf.ExtSeedSize, "ext-seed-size",
		dbConf.ExtSeedSize,
		"The additional residues to require for each seed match.")
	flag.BoolVar(&dbConf.SavePlain, "plain",
		dbConf.SavePlain,
		"When set, additional plain-text versions of files that are \n"+
			"\tnormally encoded in binary are saved with a '.plain' \n"+
			"\textension. Note that the original binary files are also saved.")
	flag.BoolVar(&dbConf.ReadOnly, "read-only",
		dbConf.ReadOnly,
		"When set, the database created will be read-only (i.e., it \n"+
			"\tcannot be appended to), but it will be smaller.")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagAppend, "append", flagAppend,
		"When set, compressed sequences will be added to existing database.")
	flag.BoolVar(&flagOverwrite, "overwrite", flagOverwrite,
		"When set, any existing database will be destroyed.")
	flag.BoolVar(&flagVerbose, "verbose", flagVerbose,
		"When set, extra output will be shown to indicate progress.")
	flag.StringVar(&flagCpuProfile, "cpuprofile", flagCpuProfile,
		"When set, a CPU profile will be written to the file specified.")
	flag.StringVar(&flagMemProfile, "memprofile", flagMemProfile,
		"When set, a memory profile will be written to the file specified.")
	flag.StringVar(&flagMemStats, "memstats", flagMemStats,
		"When set, memory statistics will be written to the file specified.")
	flag.BoolVar(&flagMemInterval, "mem-interval", flagMemInterval,
		"When set, memory profile/stats will be written at some interval.")

	flag.Usage = usage
	flag.Parse()

	runtime.GOMAXPROCS(flagGoMaxProcs)
}

func main() {
	if flag.NArg() < 2 {
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

	// If both 'append' and 'overwrite' flags are set, quit because the
	// combination doesn't make sense.
	if flagAppend && flagOverwrite {
		fatalf("Both the 'append' and 'overwrite' flags are set. It does " +
			"not make sense to set both of these flags.")
	}

	if flagVerbose {
		cablastp.Verbose = true
	}
	if flagOverwrite {
		if err := os.RemoveAll(flag.Arg(0)); err != nil {
			fatalf("Could not remove existing database '%s': %s.",
				flag.Arg(0), err)
		}
	}
	db, err := cablastp.NewWriteDB(flagAppend, dbConf, flag.Arg(0))
	if err != nil {
		fatalf("%s\n", err)
	}
	cablastp.Vprintln("")

	attachSignalHandler(db)
	pool := startCompressWorkers(db)
	orgSeqId := db.ComDB.NumSequences()
	for _, arg := range flag.Args()[1:] {
		seqChan, err := cablastp.ReadOriginalSeqs(arg, ignoredResidues)
		if err != nil {
			log.Fatal(err)
		}
		if orgSeqId == 0 {
			timer = time.Now()
		}
		for readSeq := range seqChan {
			if readSeq.Err != nil {
				log.Fatal(err)
			}
			orgSeqId = pool.compress(orgSeqId, readSeq.Seq)
			verboseOutput(db, orgSeqId)
		}
	}
	pool.done()

	cablastp.Vprintln("\n")
	cablastp.Vprintf("Wrote %s.\n", cablastp.FileCompressed)
	cablastp.Vprintf("Wrote %s.\n", cablastp.FileIndex)
	cleanup(db)
}

func cleanup(db *cablastp.DB) {
	if len(flagCpuProfile) > 0 {
		pprof.StopCPUProfile()
	}
	if len(flagMemProfile) > 0 {
		writeMemProfile(fmt.Sprintf("%s.last", flagMemProfile))
	}
	if len(flagMemStats) > 0 {
		writeMemStats(fmt.Sprintf("%s.last", flagMemStats))
	}
	if err := db.Save(); err != nil {
		fatalf("Could not save database: %s\n", err)
	}
	db.WriteClose()
}

func attachSignalHandler(db *cablastp.DB) {
	sigChan := make(chan os.Signal, 1)
	go func() {
		<-sigChan
		cleanup(db)
		os.Exit(0)
	}()
	signal.Notify(sigChan, os.Interrupt, os.Kill)
}

func verboseOutput(db *cablastp.DB, orgSeqId int) {
	if orgSeqId%interval == 0 {
		if flagVerbose {
			secElapsed := time.Since(timer).Seconds()
			seqsPerSec := float64(interval) / float64(secElapsed)

			fmt.Printf(
				"\r%d sequences compressed (%0.4f seqs/sec)",
				orgSeqId, seqsPerSec)
			timer = time.Now()
		}
		if flagMemInterval {
			if len(flagMemProfile) > 0 {
				writeMemProfile(
					fmt.Sprintf("%s.%d", flagMemProfile, orgSeqId))
			}
			if len(flagMemStats) > 0 {
				writeMemStats(
					fmt.Sprintf("%s.%d", flagMemStats, orgSeqId))
			}
		}
	}
}

func errorf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func fatalf(format string, v ...interface{}) {
	errorf(format, v...)
	os.Exit(1)
}

func writeMemProfile(name string) {
	f, err := os.Create(name)
	if err != nil {
		fatalf("%s\n", err)
	}
	pprof.WriteHeapProfile(f)
	f.Close()
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] "+
			"database-directory "+
			"fasta-file [fasta-file ...]\n",
		path.Base(os.Args[0]))
	cablastp.PrintFlagDefaults()
	os.Exit(1)
}

func writeMemStats(name string) {
	f, err := os.Create(name)
	if err != nil {
		fatalf("%s\n", err)
	}

	kb := uint64(1024)
	mb := kb * 1024

	ms := &runtime.MemStats{}
	runtime.ReadMemStats(ms)
	fmt.Fprintf(f,
		`Alloc: %d MB
TotalAlloc: %d MB
Sys: %d MB
Lookups: %d
Mallocs: %d
Frees: %d

HeapAlloc: %d MB
HeapSys: %d MB
HeapIdle: %d MB
HeapInuse: %d MB
HeapReleased: %d B
HeapObjects: %d

StackInuse: %d
StackSys: %d
MSpanInuse: %d
MSpanSys: %d
MCacheInuse: %d
MCacheSys: %d
BuckHashSys: %d

NextGC: %d
LastGC: %d
PauseTotalNs: %d s
PauseNs: %d
NumGC: %d
`,
		ms.Alloc/mb, ms.TotalAlloc/mb,
		ms.Sys/mb, ms.Lookups, ms.Mallocs,
		ms.Frees, ms.HeapAlloc/mb, ms.HeapSys/mb,
		ms.HeapIdle/mb,
		ms.HeapInuse/mb, ms.HeapReleased, ms.HeapObjects,
		ms.StackInuse, ms.StackSys, ms.MSpanInuse, ms.MSpanSys,
		ms.MCacheInuse, ms.MCacheSys, ms.BuckHashSys,
		ms.NextGC, ms.LastGC, ms.PauseTotalNs/1000000000,
		ms.PauseNs, ms.NumGC)

	f.Close()
}
