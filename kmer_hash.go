package main

import (
	"log"
)

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

// aminoValue gets the base-20 index of a particular resiude.
// aminoValue assumes that letter is a valid ASCII character in the
// inclusive range 'A' ... 'Z'.
// If the value returned by the alphaNums map is -1, aminoValue panics.
func aminoValue(letter byte) int {
	val := alphaNums[letter - 'A']
	if val == -1 {
		log.Panicf("Invalid amino acid letter: %s", b)
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
