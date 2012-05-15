package main

import (
	"log"
)

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

func hashKmer(kmer []byte) int {
	key := 0
	for i, b := range kmer {
		val := alphaNums[b - 'A']
		if val == -1 {
			log.Panicf("Invalid amino acid letter: %s", b)
		}

		key += val * pow(alphaSize, i)
	}
	return key
}
