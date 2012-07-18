package main

import "log"

// alphaSize is the number of letters in the K-mer alphabet.
const alphaSize = 20

// alphaNums is a map to assign *valid* amino acid resiudes contiunous values
// so that base-20 arithmetic can be performed on them.
// Invalid amino acid resiudes map to -1 and will produce a panic.
var alphaNums = []int{
	0,  // 'A'
	-1, // 'B'
	1,  // 'C'
	2,  // 'D'
	3,  // 'E'
	4,  // 'F'
	5,  // 'G'
	6,  // 'H'
	7,  // 'I'
	-1, // 'J'
	8,  // 'K'
	9,  // 'L'
	10, // 'M'
	11, // 'N'
	-1, // 'O'
	12, // 'P'
	13, // 'Q'
	14, // 'R'
	15, // 'S'
	16, // 'T'
	-1, // 'U'
	17, // 'V'
	18, // 'W'
	-1, // 'X'
	19, // 'Y'
	-1, // 'Z'
}

// seedLoc represents the information required to translate a seed to a slice
// of residues from the reference database. Namely, the index of the sequence
// in the reference database and the index of the residue where the seed starts
// in that sequence. The length of the seed is a constant set at
// run-time: flagSeedSize.
type seedLoc struct {
	seqInd, resInd int
}

func newSeedLoc(seqInd, resInd int) seedLoc {
	return seedLoc{
		seqInd: seqInd,
		resInd: resInd,
	}
}

func (loc seedLoc) String() string {
	return fmt.Sprintf("(%d, %d)", loc.seqInd, loc.resInd)
}

type seeds [][]seedLoc

// newSeeds creates a new table of seed location lists. The table is
// initialized with enough memory to hold lists for all possible K-mers.
// Namely, the length of seeds is equivalent to 20^(K) where 20 is the number
// of amino acids (size of alphabet) and K is equivalent to the length of
// each K-mer.
func newSeeds() seeds {
	return make(seeds, pow(alphaSize, flagSeedSize))
}

// add will create seed locations for all K-mers in refSeq and add them to
// the seeds table. Invalid K-mers are automatically skipped.
func (ss seeds) add(refSeqIndex int, refSeq *cablastp.ReferenceSeq) {
	for i := 0; i < refSeq.Len()-flagSeedSize; i++ {
		kmer := refSeq.residues[i : i+flagSeedSize]
		if !allUpperAlpha(kmer) {
			continue
		}

		kmerIndex := hashKmer(kmer)
		loc := newSeedLoc(refSeqIndex, i)

		// If no memory has been allocated for this kmer, then do so now.
		if ss[kmerIndex] == nil {
			ss[kmerIndex] = make([]seedLoc, 1)
			ss[kmerIndex][0] = loc
		} else {
			ss[kmerIndex] = append(ss[kmerIndex], loc)
		}
	}
}

// lookup returns a list of all seed locations corresponding to a particular
// K-mer.
func (ss seeds) lookup(kmer []byte) []seedLoc {
	return ss[hashKmer(kmer)]
}

// String returns a fairly crude visual representation of the seeds table.
func (ss seeds) String() string {
	strs := make([]string, 0, len(ss))
	for key, locs := range ss {
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
// If the value returned by the alphaNums map is -1, aminoValue panics.
func aminoValue(letter byte) int {
	val := alphaNums[letter-'A']
	if val == -1 {
		log.Panicf("Invalid amino acid letter: %s", letter)
	}
	return val
}

// hashKmer returns a unique hash of any 'kmer'. hashKmer assumes that 'kmer'
// satisfies 'allUpperAlpha'. hashKmer panics if there are any invalid amino
// acid residues (i.e., 'J').
//
// hashKmer satisfies this law:
// Forall a, b in [A..Z]*, hashKmer(a) == hashKmer(b) IFF a == b.
func hashKmer(kmer []byte) int {
	key := 0
	for i, b := range kmer {
		key += aminoValue(b) * pow(alphaSize, i)
	}
	return key
}

// allUpperAlpha returns true if all values in the 'kmer' slice correspond to
// values in the inclusive ASCII range 'A' ... 'Z'.
func allUpperAlpha(kmer []byte) bool {
	for _, b := range kmer {
		i := int(b - 'A')
		if i < 0 || i >= len(alphaNums) {
			return false
		}
	}
	return true
}
