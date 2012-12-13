package cablastp

import (
	"log"
	"strings"
	"sync"

	"github.com/BurntSushi/cablastp/blosum"
)

// SeedAlphaNums is a map to assign *valid* amino acid resiudes contiunous
// values so that base-N arithmetic can be performed on them. (Where
// N = SeedAlphaSize.)
// Invalid amino acid resiudes map to -1 and will produce a panic.
var (
	SeedAlphaSize        = len(blosum.Alphabet62)
	SeedAlphaNums        = make([]int, 26)
	ReverseSeedAlphaNums = make([]byte, 26)
)

// Populate SeedAlphaNums and ReverseSeedAlphaNums using the BLOSUM62
// alphabet data.
func init() {
	var amino byte

	aminoVal := 0
	for i := byte(0); i < 26; i++ {
		amino = 'A' + i
		if strings.ContainsRune(blosum.Alphabet62, rune(amino)) {
			SeedAlphaNums[i] = aminoVal
			ReverseSeedAlphaNums[aminoVal] = amino
			aminoVal += 1
		} else {
			SeedAlphaNums[i] = -1
		}
	}
}

// SeedLoc represents the information required to translate a seed to a slice
// of residues from the coarse database. Namely, the index of the sequence
// in the coarse database and the index of the residue where the seed starts
// in that sequence.
//
// Every SeedLoc also contains a pointer to the next seed location. This
// design was chosen so that each SeedLoc is independently allocated (as
// opposed to using a slice, which incurs a lot of allocation overhead
// when expanding the slice, and has the potential for pinning memory).
type SeedLoc struct {
	// Index into the coarse database sequence slice.
	SeqInd uint32

	// Index into the coarse sequence corresponding to `SeqInd`.
	ResInd uint16

	Next *SeedLoc
}

func NewSeedLoc(seqInd uint32, resInd uint16) *SeedLoc {
	return &SeedLoc{seqInd, resInd, nil}
}

// Seeds is a list of lists of seed locations. The index into the seeds table
// corresponds to a hash of particular K-mer. The list found at each
// row in the seed table corresponds to all locations in the coarse database
// in which the K-mer occurs.
type Seeds struct {
	// Table of lists of seed locations. Its length is always equivalent
	// to (SeedAlphaSize)^(SeedSize).
	Locs []*SeedLoc

	SeedSize            int
	lowComplexityWindow int // The low complexity region window size.

	// Lock is used to make Seeds.Add and Seeds.Lookup thread safe.
	lock *sync.RWMutex

	// Cache
	// (SeedAlphaSize)^0, (SeedAlphaSize)^1 ... (SeedAlphaSize)^(SeedAlphaSize).
	powers []int

	// The total number of seeds in the table.
	numSeeds int64
}

// NewSeeds creates a new table of seed location lists. The table is
// initialized with enough memory to hold lists for all possible K-mers.
// Namely, the length of seeds is equivalent to 20^(K) where 20 is the number
// of amino acids (size of alphabet) and K is equivalent to the length of
// each K-mer.
func NewSeeds(seedSize, lowComplexityWindow int) Seeds {
	powers := make([]int, seedSize+1)
	p := 1
	for i := 0; i < len(powers); i++ {
		powers[i] = p
		p *= SeedAlphaSize
	}

	locs := make([]*SeedLoc, powers[seedSize])

	return Seeds{
		Locs:                locs,
		SeedSize:            seedSize,
		lowComplexityWindow: lowComplexityWindow,
		lock:                &sync.RWMutex{},
		powers:              powers,
		numSeeds:            0,
	}
}

// NumSeeds computes the number of seeds currently in the seeds table.
// Since the seeds table is typically big, this is an expensive operation.
func (ss Seeds) NumSeeds() int64 {
	ss.lock.RLock()
	defer ss.lock.RUnlock()

	return ss.numSeeds
}

// MaybeWipe completely wipes the seeds table if the memory of the seeds table
// exceeds seedTableSizeGB (which is the number of gigabytes).
func (ss *Seeds) MaybeWipe(seedTableSizeGB float64) {
	seedLocSize := int64(16)
	maxSeedBytes := seedTableSizeGB * 1024.0 * 1024.0 * 1024.0
	maxSeeds := int64(maxSeedBytes) / seedLocSize
	if ss.NumSeeds() >= maxSeeds {
		println("Blowing away seeds table...")
		ss.lock.Lock() // acquire write lock to blow away seeds table
		defer ss.lock.Unlock()

		for i := range ss.Locs {
			ss.Locs[i] = nil
		}
		ss.numSeeds = 0
	}
}

// Add will create seed locations for all K-mers in corSeq and add them to
// the seeds table.
func (ss *Seeds) Add(coarseSeqIndex int, corSeq *CoarseSeq) {
	ss.lock.Lock()
	// Don't use defer. It comes with a performance penalty in hot spots.

	for i := 0; i < corSeq.Len()-ss.SeedSize; i++ {
		if IsLowComplexity(corSeq.Residues, i, ss.lowComplexityWindow) {
			continue
		}

		kmer := corSeq.Residues[i : i+ss.SeedSize]

		kmerIndex := ss.hashKmer(kmer)
		loc := NewSeedLoc(uint32(coarseSeqIndex), uint16(i))
		ss.numSeeds++

		if ss.Locs[kmerIndex] == nil {
			ss.Locs[kmerIndex] = loc
		} else {
			lk := ss.Locs[kmerIndex]
			for ; lk.Next != nil; lk = lk.Next {
			}
			lk.Next = loc
		}
	}

	ss.lock.Unlock()
}

// Lookup returns a list of all seed locations corresponding to a particular
// K-mer.
//
// `mem` is a pointer to a slice of seed locations, where a seed location is
// a tuple of (sequence index, residue index). `mem` is used to prevent
// unnecessary allocation. A pointer to thise slice is returned.
func (ss Seeds) Lookup(kmer []byte, mem *[][2]uint) [][2]uint {
	ss.lock.RLock()
	// Don't use defer. It comes with a performance penalty in hot spots.

	seeds := ss.Locs[ss.hashKmer(kmer)]
	if seeds == nil {
		ss.lock.RUnlock()
		return nil
	}
	*mem = (*mem)[:0]
	for seedLoc := seeds; seedLoc != nil; seedLoc = seedLoc.Next {
		*mem = append(*mem,
			[2]uint{uint(seedLoc.SeqInd), uint(seedLoc.ResInd)})
	}
	ss.lock.RUnlock()

	return *mem
}

// aminoValue gets the base-20 index of a particular resiude.
// aminoValue assumes that letter is a valid ASCII character in the
// inclusive range 'A' ... 'Z'.
// If the value returned by the SeedAlphaNums map is -1, aminoValue panics.
// (A value of -1 in this case corresponds to a residue that isn't in
// BLOSUM62.)
func aminoValue(letter byte) int {
	val := SeedAlphaNums[letter-'A']
	if val == -1 {
		log.Panicf("Invalid amino acid letter: %c", letter)
	}
	return val
}

// hashKmer returns a unique hash of any 'kmer'. hashKmer assumes that 'kmer'
// contains only valid upper case residue characters.
// hashKmer panics if there are any invalid amino acid residues.
//
// hashKmer satisfies this law:
// Forall a, b in [A..Z]*, hashKmer(a) == hashKmer(b) IFF a == b.
func (ss Seeds) hashKmer(kmer []byte) int {
	hash := 0
	lastPow := len(kmer) - 1
	for i, b := range kmer {
		hash += SeedAlphaNums[b-'A'] * ss.powers[lastPow-i]
	}
	return hash
}

// unhashKmer reverses hashKmer. (This is used in decompression.)
func (ss Seeds) unhashKmer(hash int) []byte {
	residues := make([]byte, 0)
	base := SeedAlphaSize
	for i := 0; i < ss.SeedSize; i++ {
		onesZeroed := (hash / base) * base
		digit := hash - onesZeroed
		residues = append(residues, ReverseSeedAlphaNums[digit])
		hash /= base
	}
	for i, j := 0, len(residues)-1; i < j; i, j = i+1, j-1 {
		residues[i], residues[j] = residues[j], residues[i]
	}
	return residues
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
