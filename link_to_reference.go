package cablastp

import (
	"fmt"
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
	Diff             string
	RefSeqId         int
	RefStart, RefEnd int
}

func NewLinkToReference(refSeqId, refStart, refEnd int,
	alignment [2][]byte) LinkToReference {

	return LinkToReference{
		Diff:     NewEditScript(alignment).String(),
		RefSeqId: refSeqId,
		RefStart: refStart,
		RefEnd:   refEnd,
	}
}

func NewLinkToReferenceNoDiff(refSeqId, refStart, refEnd int) LinkToReference {
	return LinkToReference{
		Diff:     "",
		RefSeqId: refSeqId,
		RefStart: refStart,
		RefEnd:   refEnd,
	}
}

func (lk LinkToReference) String() string {
	return fmt.Sprintf(
		"reference sequence id: %d, reference range: (%d, %d)\n%s",
		lk.RefSeqId, lk.RefStart, lk.RefEnd, lk.Diff)
}
