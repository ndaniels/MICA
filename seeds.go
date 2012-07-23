package cablastp

import (
	"fmt"
	"log"
	"math"
	"strings"
	"sync"
)

// SeedAlphaSize is the number of letters in the K-mer alphabet.
const SeedAlphaSize = 23

// SeedAlphaNums is a map to assign *valid* amino acid resiudes contiunous 
// values so that base-26 arithmetic can be performed on them.
// Invalid amino acid resiudes map to -1 and will produce a panic.
var SeedAlphaNums = []int{
	0,  // 'A'
	1,  // 'B'
	2,  // 'C'
	3,  // 'D'
	4,  // 'E'
	5,  // 'F'
	6,  // 'G'
	7,  // 'H'
	8,  // 'I'
	-1, // 'J'
	9,  // 'K'
	10, // 'L'
	11, // 'M'
	12, // 'N'
	-1, // 'O'
	13, // 'P'
	14, // 'Q'
	15, // 'R'
	16, // 'S'
	17, // 'T'
	-1, // 'U'
	18, // 'V'
	19, // 'W'
	20, // 'X'
	21, // 'Y'
	22, // 'Z'
}

// SeedLoc represents the information required to translate a seed to a slice
// of residues from the reference database. Namely, the index of the sequence
// in the reference database and the index of the residue where the seed starts
// in that sequence. The length of the seed is a constant set at
// run-time: flagSeedSize.
// type SeedLoc struct { 
// SeqInd, ResInd int 
// } 
type SeedLoc [2]int

func NewSeedLoc(seqInd, resInd int) SeedLoc {
	return SeedLoc{seqInd, resInd}
}

func (loc SeedLoc) String() string {
	return fmt.Sprintf("(%d, %d)", loc[0], loc[1])
}

type Seeds struct {
	Locs     [][]SeedLoc
	SeedSize int
	lock     *sync.RWMutex
}

// newSeeds creates a new table of seed location lists. The table is
// initialized with enough memory to hold lists for all possible K-mers.
// Namely, the length of seeds is equivalent to 20^(K) where 20 is the number
// of amino acids (size of alphabet) and K is equivalent to the length of
// each K-mer.
func NewSeeds(seedSize int) *Seeds {
	return &Seeds{
		Locs:     make([][]SeedLoc, pow(SeedAlphaSize, seedSize)),
		SeedSize: seedSize,
		lock:     &sync.RWMutex{},
	}
}

// add will create seed locations for all K-mers in refSeq and add them to
// the seeds table. Invalid K-mers are automatically skipped.
func (ss *Seeds) Add(refSeqIndex int, refSeq *ReferenceSeq) {
	for i := 0; i < refSeq.Len()-ss.SeedSize; i++ {
		kmer := refSeq.Residues[i : i+ss.SeedSize]
		if !KmerAllUpperAlpha(kmer) {
			continue
		}

		kmerIndex := hashKmer(kmer)
		loc := NewSeedLoc(refSeqIndex, i)

		// If no memory has been allocated for this kmer, then do so now.
		ss.lock.Lock()
		if ss.Locs[kmerIndex] == nil {
			ss.Locs[kmerIndex] = make([]SeedLoc, 1)
			ss.Locs[kmerIndex][0] = loc
		} else {
			ss.Locs[kmerIndex] = append(ss.Locs[kmerIndex], loc)
		}
		ss.lock.Unlock()
	}
}

// lookup returns a list of all seed locations corresponding to a particular
// K-mer.
func (ss *Seeds) Lookup(kmer []byte) []SeedLoc {
	ss.lock.RLock()
	defer ss.lock.RUnlock()

	seeds := ss.Locs[hashKmer(kmer)]
	cpy := make([]SeedLoc, len(seeds))
	copy(cpy, seeds)
	return cpy
}

// String returns a fairly crude visual representation of the seeds table.
func (ss *Seeds) String() string {
	ss.lock.RLock()
	defer ss.lock.RUnlock()

	strs := make([]string, 0, len(ss.Locs))
	for key, locs := range ss.Locs {
		if locs == nil {
			continue
		}

		lstrs := make([]string, len(locs))
		for i, loc := range locs {
			lstrs[i] = loc.String()
		}

		strs = append(strs,
			fmt.Sprintf("%d: %s", key, strings.Join(lstrs, " ")))
	}
	return strings.Join(strs, "\n")
}

// aminoValue gets the base-20 index of a particular resiude.
// aminoValue assumes that letter is a valid ASCII character in the
// inclusive range 'A' ... 'Z'.
// If the value returned by the SeedAlphaNums map is -1, aminoValue panics.
func aminoValue(letter byte) int {
	val := SeedAlphaNums[letter-'A']
	if val == -1 {
		log.Panicf("Invalid amino acid letter: %c", letter)
	}
	return val
}

// hashKmer returns a unique hash of any 'kmer'. hashKmer assumes that 'kmer'
// satisfies 'KmerAllUpperAlpha'. hashKmer panics if there are any invalid amino
// acid residues (i.e., 'J').
//
// hashKmer satisfies this law:
// Forall a, b in [A..Z]*, hashKmer(a) == hashKmer(b) IFF a == b.
func hashKmer(kmer []byte) int {
	key := 0
	for i, b := range kmer {
		key += aminoValue(b) * pow(SeedAlphaSize, i)
	}
	return key
}

// KmerAllUpperAlpha returns true if all values in the 'kmer' slice correspond 
// to values in the inclusive ASCII range 'A' ... 'Z'.
func KmerAllUpperAlpha(kmer []byte) bool {
	for _, b := range kmer {
		i := int(b - 'A')
		if i < 0 || i >= len(SeedAlphaNums) {
			return false
		}
	}
	return true
}

func pow(x, y int) int {
	return int(math.Pow(float64(x), float64(y)))
}
