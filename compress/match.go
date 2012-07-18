package main

import (
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
	return len(m1.orgRes) < len(m2.orgRes)
}
