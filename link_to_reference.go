package cablastp

import (
	"code.google.com/p/biogo/seq"
)

// LinkToReference represents a component of a compressed original sequence
// that allows perfect reconstruction (i.e., decompression) of the original
// sequence.
type LinkToReference struct {
	// Diff, when "applied" to the porition of the reference sequence indicated
	// by this link, will yield the original sequence corresponding to this
	// link precisely. If Diff is empty, then the subsequence of the reference
	// sequence indicated here is equivalent to the corresponding piece of
	// the original sequence.
	Diff             EditScript
	RefSeqId         int
	RefStart, RefEnd int
}

func NewLinkToReference(refSeqId, refStart, refEnd int,
	alignment seq.Alignment) *LinkToReference {

	return &LinkToReference{
		Diff:     NewEditScript(alignment),
		RefSeqId: refSeqId,
		RefStart: refStart,
		RefEnd:   refEnd,
	}
}

func NewLinkToReferenceNoDiff(refSeqId, refStart, refEnd int) *LinkToReference {
	return &LinkToReference{
		Diff:     EditScript(""),
		RefSeqId: refSeqId,
		RefStart: refStart,
		RefEnd:   refEnd,
	}
}

type EditScript string

func NewEditScript(alignment seq.Alignment) EditScript {
	return ""
}
