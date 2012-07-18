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

func NewLinkToReference(refSeqId int, refSeq *ReferenceSeq,
	alignment seq.Alignment) *LinkToReference {

	return &LinkToReference{
		Diff:     NewEditScript(alignment),
		RefSeqId: refSeqId,
		RefStart: refSeq.Offset,
		RefEnd:   refSeq.Offset + refSeq.Len(),
	}
}

func NewLinkToReferenceSimple(
	refSeqId int, refSeq *ReferenceSeq) *LinkToReference {

	return &LinkToReference{
		Diff:     EditScript(""),
		RefSeqId: refSeqId,
		RefStart: refSeq.Offset,
		RefEnd:   refSeq.Offset + refSeq.Len(),
	}
}

type EditScript string

func NewEditScript(alignment seq.Alignment) EditScript {
	return ""
}
