package cablastp

import (
	"fmt"
)

// LinkToCoarse represents a component of a compressed original sequence
// that allows perfect reconstruction (i.e., decompression) of the original
// sequence.
type LinkToCoarse struct {
	// Diff, when "applied" to the porition of the reference sequence indicated
	// by this link, will yield the original sequence corresponding to this
	// link precisely. If Diff is empty, then the subsequence of the reference
	// sequence indicated here is equivalent to the corresponding piece of
	// the original sequence.
	Diff                   string
	CoarseSeqId            uint
	CoarseStart, CoarseEnd uint16
}

func NewLinkToCoarse(coarseSeqId, coarseStart, coarseEnd uint,
	alignment [2][]byte) LinkToCoarse {

	return LinkToCoarse{
		Diff:        NewEditScript(alignment).String(),
		CoarseSeqId: coarseSeqId,
		CoarseStart: uint16(coarseStart),
		CoarseEnd:   uint16(coarseEnd),
	}
}

func NewLinkToCoarseNoDiff(
	coarseSeqId, coarseStart, coarseEnd uint) LinkToCoarse {

	return LinkToCoarse{
		Diff:        "",
		CoarseSeqId: coarseSeqId,
		CoarseStart: uint16(coarseStart),
		CoarseEnd:   uint16(coarseEnd),
	}
}

func (lk LinkToCoarse) String() string {
	return fmt.Sprintf(
		"reference sequence id: %d, reference range: (%d, %d)\n%s",
		lk.CoarseSeqId, lk.CoarseStart, lk.CoarseEnd, lk.Diff)
}
