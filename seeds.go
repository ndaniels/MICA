package cablastp

import (
	"fmt"
	"log"
	"strings"
	"sync"

	"github.com/BurntSushi/cablastp/blosum"
)

// SeedAlphaNums is a map to assign *valid* amino acid resiudes contiunous 
// values so that base-26 arithmetic can be performed on them.
// Invalid amino acid resiudes map to -1 and will produce a panic.
var (
	SeedAlphaSize        = len(blosum.Alphabet62)
	SeedAlphaNums        = make([]int, 26)
	ReverseSeedAlphaNums = make([]byte, 26)
)

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
// of residues from the reference database. Namely, the index of the sequence
// in the reference database and the index of the residue where the seed starts
// in that sequence. The length of the seed is a constant set at
// run-time: flagSeedSize.
type SeedLoc struct {
	SeqInd int
	ResInd int16
	Next   *SeedLoc
}

func NewSeedLoc(seqInd, resInd int) *SeedLoc {
	return &SeedLoc{seqInd, int16(resInd), nil}
}

func (loc SeedLoc) String() string {
	return fmt.Sprintf("(%d, %d)", loc.SeqInd, loc.ResInd)
}

type Seeds struct {
	Locs     []*SeedLoc
	SeedSize int
	lock     *sync.RWMutex
	powers   []int
}

// newSeeds creates a new table of seed location lists. The table is
// initialized with enough memory to hold lists for all possible K-mers.
// Namely, the length of seeds is equivalent to 20^(K) where 20 is the number
// of amino acids (size of alphabet) and K is equivalent to the length of
// each K-mer.
func NewSeeds(seedSize int) Seeds {
	powers := make([]int, seedSize+1)
	p := 1
	for i := 0; i < len(powers); i++ {
		powers[i] = p
		p *= SeedAlphaSize
	}

	locs := make([]*SeedLoc, powers[seedSize])

	return Seeds{
		Locs:     locs,
		SeedSize: seedSize,
		lock:     &sync.RWMutex{},
		powers:   powers,
	}
}

func (ss Seeds) Size() (int, int) {
	return 0, 0
}

// add will create seed locations for all K-mers in refSeq and add them to
// the seeds table. Invalid K-mers are automatically skipped.
func (ss Seeds) Add(coarseSeqIndex int, corSeq *CoarseSeq) {
	ss.lock.Lock()

	for i := 0; i < corSeq.Len()-ss.SeedSize; i++ {
		kmer := corSeq.Residues[i : i+ss.SeedSize]
		if !KmerAllUpperAlpha(kmer) {
			continue
		}

		kmerIndex := ss.hashKmer(kmer)
		loc := NewSeedLoc(coarseSeqIndex, i)

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

// lookup returns a list of all seed locations corresponding to a particular
// K-mer.
func (ss Seeds) Lookup(kmer []byte) [][2]int {
	ss.lock.RLock()
	seeds := ss.Locs[ss.hashKmer(kmer)]
	if seeds == nil {
		ss.lock.RUnlock()
		return nil
	}
	cpy := make([][2]int, 0, 10)
	for seedLoc := seeds; seedLoc != nil; seedLoc = seedLoc.Next {
		cpy = append(cpy, [2]int{seedLoc.SeqInd, int(seedLoc.ResInd)})
	}
	ss.lock.RUnlock()

	return cpy
}

func (seedLoc *SeedLoc) Copy() *SeedLoc {
	if seedLoc.Next == nil {
		return &SeedLoc{seedLoc.SeqInd, seedLoc.ResInd, nil}
	}
	return &SeedLoc{seedLoc.SeqInd, seedLoc.ResInd, seedLoc.Next.Copy()}
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
func (ss Seeds) hashKmer(kmer []byte) int {
	hash := 0
	lastPow := len(kmer) - 1
	for i, b := range kmer {
		hash += aminoValue(b) * ss.powers[lastPow-i]
	}
	return hash
}

// unhashKmer reverses hashKmer.
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
