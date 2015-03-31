package neutronium

import (
	"strings"
	"sync"

	"github.com/BurntSushi/cablastp/blosum"
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

// Seeds is a list of lists of seed locations. The index into the seeds table
// corresponds to a hash of particular K-mer. The list found at each
// row in the seed table corresponds to all locations in the coarse database
// in which the K-mer occurs.
type SeedTable struct {
	// Table of lists of seed locations. Its length is always equivalent
	// to (SeedAlphaSize)^(SeedSize).
	Locs []map[int]bool

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
func NewSeedTable(seedSize, lowComplexityWindow int) SeedTable {
	powers := make([]int, seedSize+1)
	p := 1
	for i := 0; i < len(powers); i++ {
		powers[i] = p
		p *= SeedAlphaSize
	}

	locs := make([]map[int]bool, powers[seedSize])

	return SeedTable{
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
func (ss SeedTable) NumSeeds() int64 {
	ss.lock.RLock()
	defer ss.lock.RUnlock()

	return ss.numSeeds
}

// MaybeWipe completely wipes the seeds table if the memory of the seeds table
// exceeds seedTableSizeGB (which is the number of gigabytes).
func (ss *SeedTable) MaybeWipe(seedTableSizeGB float64) {
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

func (ss *SeedTable) Add(coarseSeqIndex int, corSeq *CoarseSeq) {
	ss.lock.Lock()
	// Don't use defer. It comes with a performance penalty in hot spots.

	for i := 0; i < corSeq.Len()-ss.SeedSize; i++ {
		if IsLowComplexity(corSeq.Residues, i, ss.lowComplexityWindow) {
			continue
		}

		kmer := corSeq.Residues[i : i+ss.SeedSize]

		kmerIndex := ss.hashKmer(kmer)

		ss.numSeeds++

		if ss.Locs[kmerIndex] == nil {

			ss.Locs[kmerIndex] = make(map[int]bool)
		}

		ss.Locs[kmerIndex][coarseSeqIndex] = true
	}

	ss.lock.Unlock()
}

func (ss SeedTable) Lookup(kmer []byte, coarseSeqIndex int) bool {
	ss.lock.RLock()
	// Don't use defer. It comes with a performance penalty in hot spots.

	kmerIndex := ss.hashKmer(kmer)
	match := ss.Locs[kmerIndex][coarseSeqIndex]

	ss.lock.RUnlock()

	return match
}

// hashKmer returns a unique hash of any 'kmer'. hashKmer assumes that 'kmer'
// contains only valid upper case residue characters.
// hashKmer panics if there are any invalid amino acid residues.
//
// hashKmer satisfies this law:
// Forall a, b in [A..Z]*, hashKmer(a) == hashKmer(b) IFF a == b.
func (ss SeedTable) hashKmer(kmer []byte) int {
	hash := 0
	lastPow := len(kmer) - 1
	for i, b := range kmer {
		hash += SeedAlphaNums[b-'A'] * ss.powers[lastPow-i]
	}
	return hash
}
