package cablastp

import (
	"code.google.com/p/biogo/seq"
)

// CompressedSeq corresponds to the components of a compressed sequence.
type CompressedSeq struct {
	// Name is an uncompressed string from the original FASTA header.
	Name string

	// Links is an ordered lists of links to portions of the reference
	// database. When all links are followed, the concatenation of each
	// sequence correspond to each link equals the entire original sequence.
	Links []*LinkToReference
}

// LinkToReference represents a component of a compressed original sequence
// that allows perfect reconstruction (i.e., decompression) of the original
// sequence.
type LinkToReference struct {
	// Diff, when "applied" to the porition of the reference sequence indicated
	// by this link, will yield the original sequence corresponding to this
	// link precisely. If Diff is empty, then the subsequence of the reference
	// sequence indicated here is equivalent to the corresponding piece of
	// the original sequence.
	Diff EditScript
	RefSeqId int
	RefStart, RefEnd int
}

func NewLinkToReference(refSeqId int, refSeq *ReferenceSeq,
	alignment seq.Alignment) *LinkToReference {

	return &LinkToReference{
		Diff: NewEditScript(alignment),
		RefSeqId: refSeqId,
		RefStart: refSeq.Offset,
		RefEnd: refSeq.Offset + refSeq.Len(),
	}
}

type EditScript string

func NewEditScript(alignment seq.Alignment) EditScript {
	return ""
}

