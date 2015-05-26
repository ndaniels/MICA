package cablastp

import "fmt"

// LinkToCompressed represents a link from a reference sequence to a
// compressed original sequence. It serves as a bridge from a BLAST hit in
// the coarse database to the corresponding original sequence that is
// redundant to the specified residue range in the reference sequence.
type LinkToCompressed struct {
	OrgSeqId               uint32
	CoarseStart, CoarseEnd uint16
	Next                   *LinkToCompressed
}

func NewLinkToCompressed(
	orgSeqId uint32, coarseStart, coarseEnd uint16) *LinkToCompressed {

	return &LinkToCompressed{
		OrgSeqId:    orgSeqId,
		CoarseStart: coarseStart,
		CoarseEnd:   coarseEnd,
		Next:        nil,
	}
}

func (lk LinkToCompressed) String() string {
	return fmt.Sprintf("original sequence id: %d, coarse range: (%d, %d)",
		lk.OrgSeqId, lk.CoarseStart, lk.CoarseEnd)
}
