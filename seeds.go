package cablastp

import (
	"fmt"
	"log"
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
func (ss Seeds) Add(refSeqIndex int, refSeq *ReferenceSeq) {
	ss.lock.Lock()

	for i := 0; i < refSeq.Len()-ss.SeedSize; i++ {
		kmer := refSeq.Residues[i : i+ss.SeedSize]
		if !KmerAllUpperAlpha(kmer) {
			continue
		}

		kmerIndex := ss.hashKmer(kmer)
		loc := NewSeedLoc(refSeqIndex, i)

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
	key := 0
	for i, b := range kmer {
		key += aminoValue(b) * ss.powers[i]
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
