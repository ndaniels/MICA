package mica

import (
	"bytes"
	"fmt"
	"github.com/ndaniels/mica/blosum"
	"os"
  // "runtime"
	"sync"
)

// redCompressPool represents a pool of workers where each worker is responsible
// for compressing a single sequence at a time.
type redCompressPool struct {
	db     *DB
	jobs   chan redCompressJob
	wg     *sync.WaitGroup
	closed bool
}

// redCompressJob values are messages sent to the pool of workers when a new
// sequence should be compressed.
type redCompressJob struct {
	redSeqId int
	redSeq   *ReducedSeq
}

// startCompressWorkers initializes a pool of compression workers.
//
// The compressPool returned can be used to compress sequences concurrently.
func StartCompressReducedWorkers(db *DB) redCompressPool {
	wg := &sync.WaitGroup{}
	jobs := make(chan redCompressJob, 200)
	pool := redCompressPool{
		db:     db,
		jobs:   jobs,
		wg:     wg,
		closed: false,
	}
  // for queries, due to inherent similarities and proximity,
  // this should only use one thread.
	for i := 0; i < 1; i++ {
		wg.Add(1)
		go pool.worker()
	}
	return pool
}

// When the program ends (either by SIGTERM or when all of the input sequences
// are compressed), 'cleanup' is executed. It writes all CPU/memory profiles
// if they're enabled, waits for the compression workers to finish, saves
// the database to disk and closes all file handles.
func CleanupDB(db *DB, pool *redCompressPool) {

	pool.done()
	if err := db.Save(); err != nil {
		fatalf("Could not save database: %s\n", err)
	}
	db.WriteClose()
}

// compress will construct a redCompressJob and send it to the worker pool.
//
// compress returns the next original sequence id to be used.
func (pool redCompressPool) CompressReduced(id int, seq *ReducedSeq) int {
	pool.jobs <- redCompressJob{
		redSeqId: id,
		redSeq:   seq,
	}
	return id + 1
}

// worker is meant to be run as a goroutine. It allocates a goroutine-specific
// memory arena (to prevent allocation in hot spots like alignment and
// seed lookup), and sends the compressed sequences to the compressed
// database for writing.
func (pool redCompressPool) worker() {
	mem := newMemory()
	for job := range pool.jobs {
		comSeq := CompressReduced(pool.db, job.redSeqId, job.redSeq, mem)
		pool.db.ComDB.Write(comSeq)
	}
	pool.wg.Done()
}

// done 'joins' the worker goroutines. (Blocks until all workers are finished
// compressing sequences.)
func (pool *redCompressPool) done() {
	if pool.closed {
		return
	}
	pool.closed = true
	close(pool.jobs)
	pool.wg.Wait()
}

// TODO: need to pass around a pair of nativeSeq and reducedSeq
// What gets saved to compressed db is nativeSeq
// What gets saved to coarse fasta is reducedSeq

// CompressReduced will convert a 4-character alphabet sequence into a compressed
// sequence.
// The process involves finding commonality in the original sequence with
// other sequences in the coarse database, and linking those common
// sub-sequences to sub-sequences in the coarse database.
//
// N.B. `mem` is used in alignment and seed lookups to prevent allocation.
// Think of them as goroutine-specific memory arenas.
func CompressReduced(db *DB, redSeqId int,
	redSeq *ReducedSeq, mem *memory) CompressedSeq {

	// cseqExt and oseqExt will contain `extSeedSize` residues after the end
	// of any particular seed in coarse and original sequences, respectively.
	// If the residues are not equivalent, that particular seed is skipped.
	var cseqExt, oseqExt []byte

	// Start the creation of a compressed sequence.
	cseq := NewCompressedSeq(redSeqId, redSeq.Name)

	// Convenient aliases.
	coarsedb := db.CoarseDB
	mapSeedSize := db.MapSeedSize
	extSeedSize := db.ExtSeedSize
	olen := redSeq.Len()

	// Keep track of two pointers. 'current' refers to the residue index in the
	// original sequence that extension is currently originating from.
	// 'lastMatch' refers to the residue index of the *end* of the last match
	// with a coarse sequence in the compressed database.
	lastMatch, current := 0, 0

	// Iterate through the original sequence a 'kmer' at a time.
	skipSize := 4
	limit := olen - mapSeedSize - extSeedSize - skipSize
	for current = 0; current <= limit; current += skipSize {
		kmer := redSeq.Residues[current : current+mapSeedSize]
		// skip wildcard-containing kmers
		if bytes.IndexByte(kmer, 'N') > -1 {
			continue
		}

		seeds := coarsedb.Seeds.Lookup(kmer, &mem.seeds)

		// Before trying to extend this with seeds, check to see if there is
		// a low complexity region within `db.MinMatchLen` residues from
		// `current`. If there is, skip ahead to the end of it.
		// note that we are checking the original, unreduced sequence here!
		if db.LowComplexity > 0 {
			skip := skipLowComplexity(
				redSeq.Residues[current:], db.MinMatchLen, db.LowComplexity)
			if skip > 0 {
				current += skip
				continue
			}
		}

		// Each seed location corresponding to the current K-mer must be
		// used to attempt to extend a match.
		for _, seedLoc := range seeds {
			corSeqId := int(seedLoc[0])
			corResInd := int(seedLoc[1])
			corSeq := coarsedb.CoarseSeqGet(uint(corSeqId))

			// If the seed extension extends beyond the end of the coarse
			// sequence pointed to by seedLoc, then move along.
			extCorStart := corResInd + mapSeedSize
			extOrgStart := current + mapSeedSize
			if extCorStart+extSeedSize >= corSeq.Len() {
				continue
			}

			// If the seed extensions in each sequence are not equivalent,
			// skip this seedLoc.
			cseqExt = corSeq.Residues[extCorStart : extCorStart+extSeedSize]
			oseqExt = redSeq.Residues[extOrgStart : extOrgStart+extSeedSize]
			if !bytes.Equal(cseqExt, oseqExt) {
				continue
			}

			// The "match" between coarse and original (reduced) sequence will
			// occur somewhere between the the residue index of the seed and
			// the end of the sequence for the coarse sequence, and the
			// position of the "current" pointer and the end of the sequence
			// for the original sequence.
			// TODO maybe the simplest way to match backwards is to reverse the
			// coarse and reduced sequences FROM the last match TO the current,
			// and call extendMatch on those as well. Reverse the result, and
			// prepend it to corMatch and redMatch here.
			corMatch, redMatch := extendMatch(
				corSeq.Residues[corResInd:], redSeq.Residues[current:],
				db.GappedWindowSize, db.UngappedWindowSize,
				db.MatchKmerSize, db.ExtSeqIdThreshold,
				mem)
       

			// If the part of the original (reduced) sequence does not exceed the
			// minimum match length, then we don't accept the match and move
			// on to the next one.
			if len(redMatch) < db.MinMatchLen {
				continue
			}



			alignment := nwAlign(corMatch, redMatch, mem)
			id := SeqIdentity(alignment[0], alignment[1])
			if id < db.MatchSeqIdThreshold {
				continue
			}


			// If we end up extending a match because we're close to
			// some boundary (either a sequence or a match boundary), then
			// we need to perform another alignment.
			changed := false

			// If we're close to the end of the original sequence, extend
			// the match to the end.
			if len(redMatch)+db.MatchExtend >= redSeq.Len()-int(current) {
				redMatch = redSeq.Residues[current:]
				changed = true
			}

			// And if we're close to the end of the last match, extend this
			// match backwards.
			if current-lastMatch <= db.MatchExtend {
				end := current + len(redMatch)
				redMatch = redSeq.Residues[lastMatch:end]
				current = lastMatch
				changed = true
			}

			// If we've extended our match, we need another alignment.
			if changed {
				alignment = nwAlign(corMatch, redMatch, mem)
			}

			// Otherwise, we accept the first valid match and move on to the
			// next kmer after the match ends.
			corStart := corResInd
			corEnd := corStart + len(corMatch)
			orgStart := current
			orgEnd := orgStart + len(redMatch)

			// If there are residues between the end of the last match
			// and the start of this match, then that means no good match
			// could be found for those residues. Thus, they are added to
			// the coarse database. (A pathological LinkToCoarse is
			// created with an empty diff script that points to the added
			// region in the coarse database in its entirety.)
			if orgStart-lastMatch > 0 {
				redSub := redSeq.NewSubSequence(
					uint(lastMatch), uint(current))
				addReducedWithoutMatch(&cseq, coarsedb, redSeqId, redSub)
			}

			// For the given match, add a LinkToCoarse to the portion of
			// the coarse sequence matched. This serves as a component
			// of a compressed original sequence. Also, add a
			// LinkToCompressed to the coarse sequence matched. This
			// serves as a bridge to expand coarse sequences into their
			// original sequences.
			cseq.Add(NewLinkToCoarse(
				uint(corSeqId), uint(corStart), uint(corEnd), alignment))
			corSeq.AddLink(NewLinkToCompressed(
				uint32(redSeqId), uint16(corStart), uint16(corEnd)))

			// Skip the current pointer ahead to the end of this match.
			// Update the lastMatch pointer to point at the end of this
			// match.
			lastMatch = orgEnd
			current = orgEnd - 1

			// Don't process any more seedLocs for this K-mer once we've
			// found a match.
			break
		}
	}

	// If there are any leftover residues, then no good match for them
	// could be found. Therefore, add them to the coarse database and
	// create the appropriate links.
	if redSeq.Len()-lastMatch > 0 {
		redSub := redSeq.NewSubSequence(uint(lastMatch), uint(redSeq.Len()))
		addReducedWithoutMatch(&cseq, coarsedb, redSeqId, redSub)
	}

	return cseq
}

func reverse(a []byte) []byte {
	l := len(a)
	result := make([]byte, l)
	for i, v := range a {
		result[l-i-1] = v
	}
	return result
}

// extendMatch uses a combination of ungapped and gapped extension to find
// quality candidates for compression.
func extendMatch(corRes, orgRes []byte,
	gappedWindowSize, ungappedWindowSize, kmerSize, idThreshold int,
	mem *memory) (corMatchRes, orgMatchRes []byte) {

	// Starting at seedLoc.resInd and current (from 'compress'), corMatchLen
	// and orgMatchLen correspond to the length of the match in each of
	// the coarse and the original sequence, respectively.
	// At the end of the loop, the slices [seedLoc.resInd:corMatchLen]
	// and [current:orgMatchLen] will correspond to the match. (Again, this
	// is in the context of the inner loop in 'compress'. For this particular
	// function, corMatchLen and orgMatch start at 0, so that the matches
	// eventually returned correspond to the [:corMatchLen] and [:orgMatchLen]
	// slices.)
	corMatchLen, orgMatchLen := 0, 0
	for {
		// If the match has consumed either of the coarse or original
		// sequence, then we must quit with what we have.
		if corMatchLen == len(corRes) || orgMatchLen == len(orgRes) {
			break
		}

		// Ungapped extension returns an integer corresponding to the
		// number of residues that the match was extended by.
		matchLen := alignUngapped(
			corRes[corMatchLen:], orgRes[orgMatchLen:],
			ungappedWindowSize, kmerSize, idThreshold)

		// Since ungapped extension increases the coarse and
		// original sequence match portions equivalently, add the
		// match length to both.
		corMatchLen += matchLen
		orgMatchLen += matchLen

		// Gapped extension returns an alignment corresponding to the
		// window starting after the previous ungapped extension
		// ended plus the gapped window size. (It is bounded by the
		// length of each sequence.)
		alignment := nwAlign(
			corRes[corMatchLen:min(len(corRes), corMatchLen+gappedWindowSize)],
			orgRes[orgMatchLen:min(len(orgRes), orgMatchLen+gappedWindowSize)],
			mem)

		// If the alignment has a sequence identity below the
		// threshold, then gapped extension has failed. We therefore
		// quit and are forced to be satisfied with whatever
		// corMatchLen and orgMatchLen are set to.
		id := SeqIdentity(alignment[0], alignment[1])

		if id < idThreshold {
			break
		}

		// We live to die another day.
		// We need to add to the corMatch{Pos,Len} and orgMatch{Pos,Len}
		// just like we did for ungapped extension. However, an
		// alignment can correspond to two different sized subsequences
		// of the coarse and original sequence. Therefore, only
		// increase each by the corresponding sizes from the
		// alignment.
		corMatchLen += alignLen(alignment[0])
		orgMatchLen += alignLen(alignment[1])
	}

	return corRes[:corMatchLen], orgRes[:orgMatchLen]
}

// skipLowComplexity looks for a low complexity region starting at the
// beginning of `seq` and up to `windowSize`. If one is found, `x` is returned
// where `x` corresponds to the position of the first residue after
// the low complexity region has ended. If a low complexity region isn't
// found, `0` is returned.
//
// N.B. regionSize is the number of contiguous positions in the sequence
// that must contain the same residue in order to qualify as a low complexity
// region.
func skipLowComplexity(seq []byte, windowSize, regionSize int) int {
	upto := min(len(seq), windowSize+regionSize)
	last, repeats, i, found := byte(0), 1, 0, false
	for i = 0; i < upto; i++ {
		if seq[i] == last {
			repeats++
			if repeats >= regionSize {
				found = true
				break
			}
			continue
		}

		// The last residue isn't the same as this residue, so reset.
		last = seq[i]
		repeats = 1
	}
	if !found { // no low complexity region was found.
		return 0
	}

	// We're in a low complexity region. Consume as many residues equal
	// to `last` as possible.
	//
	// N.B. `i` is already set to where we left off in the last loop.
	for ; i < len(seq); i++ {
		if seq[i] != last { // end of low complexity region
			break
		}
	}
	return i
}

// addWithoutMatch adds a portion of an original sequence that could not be
// matched to anything in the coarse database to the coarse database.
// A LinkToCompressed is created and automatically added to the new coarse
// sequence.
//
// An appropriate link is also added to the given compressed sequence.
// TODO we need to juggle the reduced and original seq portions
func addReducedWithoutMatch(cseq *CompressedSeq,
	coarsedb *CoarseDB, redSeqId int, redSub *ReducedSeq) {

	// Explicitly copy residues to avoid pinning memory.
	redSubCpy := make([]byte, len(redSub.Residues))
	copy(redSubCpy, redSub.Residues)

	corSeqId, corSeq := coarsedb.Add(redSubCpy)
	corSeq.AddLink(
		NewLinkToCompressed(uint32(redSeqId), 0, uint16(len(redSubCpy))))

	cseq.Add(
		NewLinkToCoarseNoDiff(uint(corSeqId), 0, uint(len(redSubCpy))))
}

// func min(a, b int) int {
// 	if a < b {
// 		return a
// 	}
// 	return b
// }

// func max(a, b int) int {
// 	if a > b {
// 		return a
// 	}
// 	return b
// }

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// memory.go in cablastp-compress
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

const (
	memSeqSize       = 10000
	dynamicTableSize = memSeqSize * memSeqSize
	numSeeds         = 100
)

// memory is a goroutine-specific memory arena, specifically used in each
// compression worker goroutine. Its purpose is to reduce the amount of
// memory allocation in hot-spots: sequence alignment and seed lookup.
type memory struct {
	table    []int
	ref, org []byte
	seeds    [][2]uint
}

func newMemory() *memory {
	return &memory{
		table: make([]int, memSeqSize*memSeqSize),
		ref:   make([]byte, memSeqSize),
		org:   make([]byte, memSeqSize),
		seeds: make([][2]uint, 0, numSeeds),
	}
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// align.go in cablastp-compress
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// alignLen computes the length of a sequence in an alignment.
// (i.e., the number of residues that aren't "-".)
func alignLen(seq []byte) (length int) {
	for _, res := range seq {
		if res != '-' {
			length++
		}
	}
	return
}

// alignUngapped takes a coarse and an original sub-sequence and returns a
// length corresponding to the number of amino acids scanned by greedily
// consuming successive K-mer matches in N-mer windows.
//
// The algorithm works by attempting to find *exact* K-mer matches between the
// sequences in N-mer windows. If N residues are scanned and no K-mer match
// is found, the the current value of length is returned (which may be 0).
// If a K-mer match is found, the current value of length is set to the total
// number of amino acid residues scanned, and a search for the next K-mer match
// for the next N-mer window is started.
func alignUngapped(rseq []byte, oseq []byte,
	windowSize, kmerSize, idThreshold int) int {
    
	length, scanned, successive := 0, 0, 0
	tryNextWindow := true
	for tryNextWindow {
		tryNextWindow = false
		for i := 0; i < windowSize; i++ {
			// If we've scanned all residues in one of the sub-sequences, then
			// there is nothing left to do for ungapped extension. Therefore,
			// quit and return the number of residues scanned up until the
			// *last* match.
			if scanned >= len(rseq) || scanned >= len(oseq) {
				break
			}

			if rseq[scanned] == oseq[scanned] {
				successive++
			} else {
				successive = 0
			}

			scanned++
			if successive == kmerSize {
				// Get the residues between matches: i.e., after the last
				// match to the start of this match. But only if there is at
				// least one residue in that range.
				if (scanned-kmerSize)-length > 0 {
					id := SeqIdentity(
						rseq[length:scanned-kmerSize],
						oseq[length:scanned-kmerSize])

					// If the identity is less than the threshold, then this
					// K-mer match is no good. But keep trying until the window
					// is closed. (We "keep trying" by decrementing successive
					// matches by 1.)
					if id < idThreshold {
						successive--
						continue
					}
				}

				// If we're here, then we've found a valid match. Update the
				// length to indicate the number of residues scanned and make
				// sure we try the next Ungapped window.
				length = scanned
				successive = 0
				tryNextWindow = true
				break
			}
		}
	}
	return length
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// nw.go in cablastp-compress
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

var (
	resTrans [256]int
)

// Initialize the alignment lookup table. (i.e., translate ASCII residue
// characters to BLOSUM62 matrix indices.)
func init() {
	for i := 0; i < len(blosum.Alphabet62); i++ {
		resTrans[blosum.Alphabet62[i]] = i
	}
}

// nwAlign performs Needleman-Wunsch sequence alignment.
//
// This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the
// biogo package: code.google.com/p/biogo.
//
// It's mostly copied from its original form, but it is optimized specifically
// for  to limit allocations and to absolve the need for the
// biogo/seq.Seq type. There are several additional optimizations to limit
// functional calls and pointer deferences.
//
// Perhaps the biggest optimization; however, is constraining dynamic
// programming to only allow a limited number of gaps proportion to the
// length of the large of rseq and oseq.
func nwAlign(rseq, oseq []byte, mem *memory) [2][]byte {
	gap := len(blosum.Matrix62) - 1
	r, c := len(rseq)+1, len(oseq)+1
	off := 0

	constrained := true
	constraint := max(r, c) / 4
	if r <= 11 || c <= 11 {
		constrained = false
	}

	var table []int
	var j int
	if r*c > dynamicTableSize {
		table = make([]int, r*c)
	} else {
		table = mem.table[:r*c]
		for i := range table {
			table[i] = 0
		}
	}

	var sdiag, sup, sleft, rVal, oVal int
	gapChar := byte('-')
	matrix := blosum.Matrix62

	var i2, i3 int
	for i := 1; i < r; i++ {
		i2 = (i - 1) * c
		i3 = i * c
		for j = 1; j < c; j++ {
			if constrained && ((i-j) > constraint || (j-i) > constraint) {
				continue
			}
			rVal, oVal = resTrans[rseq[i-1]], resTrans[oseq[j-1]]

			off = i2 + (j - 1)
			sdiag = table[off] + matrix[rVal][oVal]
			sup = table[off+1] + matrix[rVal][gap]
			sleft = table[off+c] + matrix[gap][oVal]
			switch {
			case sdiag > sup && sdiag > sleft:
				table[i3+j] = sdiag
			case sup > sleft:
				table[i3+j] = sup
			default:
				table[i3+j] = sleft
			}
		}
	}

	refAln, orgAln := mem.ref[:0], mem.org[:0]

	i, j := r-1, c-1
	for i > 0 && j > 0 {
		rVal, oVal = resTrans[rseq[i-1]], resTrans[oseq[j-1]]

		sdiag = table[(i-1)*c+(j-1)] + matrix[rVal][oVal]
		sup = table[(i-1)*c+j] + matrix[gap][oVal]
		sleft = table[i*c+(j-1)] + matrix[rVal][gap]
		switch {
		case sdiag > sup && sdiag > sleft:
			i--
			j--
			refAln = append(refAln, rseq[i])
			orgAln = append(orgAln, oseq[j])
		case sup > sleft:
			i--
			refAln = append(refAln, rseq[i])
			orgAln = append(orgAln, gapChar)
		default:
			j--
			refAln = append(refAln, gapChar)
			orgAln = append(orgAln, oseq[j])
		}
	}

	for ; i > 0; i-- {
		refAln = append(refAln, rseq[i-1])
		orgAln = append(orgAln, gapChar)
	}
	for ; j > 0; j-- {
		refAln = append(refAln, gapChar)
		orgAln = append(orgAln, oseq[j-1])
	}

	if len(refAln) == len(orgAln) {
		for i, j := 0, len(refAln)-1; i < j; i, j = i+1, j-1 {
			refAln[i], refAln[j] = refAln[j], refAln[i]
			orgAln[i], orgAln[j] = orgAln[j], orgAln[i]
		}
	} else {
		for i, j := 0, len(refAln)-1; i < j; i, j = i+1, j-1 {
			refAln[i], refAln[j] = refAln[j], refAln[i]
		}
		for i, j := 0, len(orgAln)-1; i < j; i, j = i+1, j-1 {
			orgAln[i], orgAln[j] = orgAln[j], orgAln[i]
		}
	}

	return [2][]byte{refAln, orgAln}
}

func fatalf(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
	os.Exit(1)
}