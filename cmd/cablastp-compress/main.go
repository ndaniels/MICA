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

// makeblastdb -dbtype prot -input_type fasta

const interval = 1000

var (
	// Used to compute the number of sequences compressed per second.
	timer time.Time

	// Any residue in `ignoredResidues` will be replaced with an X.
	// These should correspond to the residues NOT in blosum.Alphabet62.
	ignoredResidues = []byte{'J', 'O', 'U'}

	// A default configuration.
	dbConf = cablastp.DefaultDBConf

	// Flags that affect the higher level operation of compression.
	// Flags that control algorithmic parameters are stored in `dbConf`.
	flagGoMaxProcs  = runtime.NumCPU()
	flagAppend      = false
	flagOverwrite   = false
	flagQuiet       = false
	flagMaxSeedsGB  = 8.0
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
	flag.IntVar(&dbConf.LowComplexity, "low-complexity",
		dbConf.LowComplexity,
		"The window size used to detect regions of low complexity.\n"+
			"\tLow complexity regions are repetitions of a single amino\n"+
			"\tacid residue. Low complexity regions are skipped when\n"+
			"\ttrying to extend a match.")
	flag.IntVar(&dbConf.SeedLowComplexity, "seed-low-complexity",
		dbConf.SeedLowComplexity,
		"The seed window size used to detect regions of low complexity.\n"+
			"\tLow complexity regions are repetitions of a single amino\n"+
			"\tacid residue. Low complexity regions matching this window\n"+
			"\tsize are not included in the seeds table.")
	flag.BoolVar(&dbConf.SavePlain, "plain",
		dbConf.SavePlain,
		"When set, additional plain-text versions of files that are \n"+
			"\tnormally encoded in binary are saved with a '.plain' \n"+
			"\textension. Note that the original binary files are also saved.")
	flag.BoolVar(&dbConf.ReadOnly, "read-only",
		dbConf.ReadOnly,
		"When set, the database created will be read-only (i.e., it \n"+
			"\tcannot be appended to), but it will be smaller.")
	flag.StringVar(&dbConf.BlastMakeBlastDB, "makeblastdb",
		dbConf.BlastMakeBlastDB,
		"The location of the 'makeblastdb' executable.")

	flag.IntVar(&flagGoMaxProcs, "p", flagGoMaxProcs,
		"The maximum number of CPUs that can be executing simultaneously.")
	flag.BoolVar(&flagAppend, "append", flagAppend,
		"When set, compressed sequences will be added to existing database.\n"+
			"\tThe parameters used to create the initial database are\n"+
			"\tautomatically used by default. They can still be overriden\n"+
			"\ton the command line.")
	flag.BoolVar(&flagOverwrite, "overwrite", flagOverwrite,
		"When set, any existing database will be destroyed.")
	flag.BoolVar(&flagQuiet, "quiet", flagQuiet,
		"When set, the only outputs will be errors echoed to stderr.")
	flag.Float64Var(&flagMaxSeedsGB, "max-seeds", flagMaxSeedsGB,
		"When set, the in memory seeds table will be completely erased\n"+
			"\twhen the memory used by seeds exceeds the specified number,\n"+
			"\tin gigabytes.\n"+
			"\tEach seed corresponds to 16 bytes of memory.\n"+
			"\tSetting to zero disables this behavior.")
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

	// If both 'append' and 'overwrite' flags are set, quit because the
	// combination doesn't make sense.
	if flagAppend && flagOverwrite {
		fatalf("Both the 'append' and 'overwrite' flags are set. It does " +
			"not make sense to set both of these flags.")
	}

	// If the quiet flag isn't set, enable verbose output.
	if !flagQuiet {
		cablastp.Verbose = true
	}

	// If the overwrite flag is set, remove whatever directory that may
	// already be there.
	if flagOverwrite {
		if err := os.RemoveAll(flag.Arg(0)); err != nil {
			fatalf("Could not remove existing database '%s': %s.",
				flag.Arg(0), err)
		}
	}

	// Create a new database for writing. If we're appending, we load
	// the coarse database into memory, and setup the database for writing.
	db, err := cablastp.NewWriteDB(flagAppend, dbConf, flag.Arg(0))
	if err != nil {
		fatalf("%s\n", err)
	}
	cablastp.Vprintln("")

	pool := startCompressWorkers(db)
	orgSeqId := db.ComDB.NumSequences()
	mainQuit := make(chan struct{}, 0)

	// If the process is killed, try to clean up elegantly.
	// The idea is to preserve the integrity of the database.
	attachSignalHandler(db, mainQuit, &pool)

	// Start the CPU profile after all of the data has been read.
	if len(flagCpuProfile) > 0 {
		f, err := os.Create(flagCpuProfile)
		if err != nil {
			fatalf("%s\n", err)
		}
		pprof.StartCPUProfile(f)
	}
	for _, arg := range flag.Args()[1:] {
		seqChan, err := cablastp.ReadOriginalSeqs(arg, ignoredResidues)
		if err != nil {
			log.Fatal(err)
		}
		if orgSeqId == 0 {
			timer = time.Now()
		}
		for readSeq := range seqChan {
			// Do a non-blocking receive to see if main needs to quit.
			select {
			case <-mainQuit:
				<-mainQuit // wait for cleanup to finish before exiting main.
				return
			default:
			}

			if readSeq.Err != nil {
				log.Fatal(err)
			}
			dbConf.BlastDBSize += uint64(readSeq.Seq.Len())
			orgSeqId = pool.compress(orgSeqId, readSeq.Seq)
			verboseOutput(db, orgSeqId)
			if flagMaxSeedsGB > 0 && orgSeqId%10000 == 0 {
				db.CoarseDB.Seeds.MaybeWipe(flagMaxSeedsGB)
			}
		}
	}
	cablastp.Vprintln("\n")
	cablastp.Vprintf("Wrote %s.\n", cablastp.FileCompressed)
	cablastp.Vprintf("Wrote %s.\n", cablastp.FileIndex)

	cleanup(db, &pool)
}

// When the program ends (either by SIGTERM or when all of the input sequences
// are compressed), 'cleanup' is executed. It writes all CPU/memory profiles
// if they're enabled, waits for the compression workers to finish, saves
// the database to disk and closes all file handles.
func cleanup(db *cablastp.DB, pool *compressPool) {
	if len(flagCpuProfile) > 0 {
		pprof.StopCPUProfile()
	}
	if len(flagMemProfile) > 0 {
		writeMemProfile(fmt.Sprintf("%s.last", flagMemProfile))
	}
	if len(flagMemStats) > 0 {
		writeMemStats(fmt.Sprintf("%s.last", flagMemStats))
	}
	pool.done()
	if err := db.Save(); err != nil {
		fatalf("Could not save database: %s\n", err)
	}
	db.WriteClose()
}

// Runs a goroutine to listen for SIGTERM and SIGKILL.
func attachSignalHandler(db *cablastp.DB, mainQuit chan struct{},
	pool *compressPool) {

	sigChan := make(chan os.Signal, 1)
	go func() {
		<-sigChan
		mainQuit <- struct{}{}
		cleanup(db, pool)
		mainQuit <- struct{}{}
		os.Exit(0)
	}()
	signal.Notify(sigChan, os.Interrupt, os.Kill)
}

// The output generated after each sequence is compressed (or more precisely,
// after some interval of sequences has been compressed).
func verboseOutput(db *cablastp.DB, orgSeqId int) {
	if orgSeqId%interval == 0 {
		if !flagQuiet {
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

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
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

func usage() {
	fmt.Fprintf(os.Stderr,
		"\nUsage: %s [flags] "+
			"database-directory "+
			"fasta-file [fasta-file ...]\n",
		path.Base(os.Args[0]))
	cablastp.PrintFlagDefaults()
	os.Exit(1)
}

// A nasty function to format the runtime.MemStats struct for human
// consumption.
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
