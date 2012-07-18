package main

import (
	"github.com/BurntSushi/cablastp"

	"code.google.com/p/biogo/seq"
)

type match struct {
	refSeqId         int
	refSeq           *cablastp.ReferenceSeq
	refStart, refEnd int
	orgStart, orgEnd int
	alignment        seq.Alignment
}

func (m1 match) Less(m2 match) bool {
	return (m1.orgEnd - m1.orgStart) < (m2.orgEnd - m2.orgStart)
}

// bestMatch searches a list of matches and returns the "best" one according
// to the definition of Less.
func bestMatch(matches []match) match {
	bestMatch := matches[0]
	for _, match := range matches[1:] {
		if bestMatch.Less(match) {
			bestMatch = match
		}
	}
	return bestMatch
}
