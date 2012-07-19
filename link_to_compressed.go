package cablastp

import "fmt"

// LinkToCompressed represents a link from a reference sequence to a
// compressed original sequence. It serves as a bridge from a BLAST hit in
// the coarse database to the corresponding original sequence that is
// redundant to the specified residue range in the reference sequence.
type LinkToCompressed struct {
	OrgSeqId         int
	RefStart, RefEnd int
}

func NewLinkToCompressed(orgSeqId, refStart, refEnd int) *LinkToCompressed {
	return &LinkToCompressed{
		OrgSeqId: orgSeqId,
		RefStart: refStart,
		RefEnd:   refEnd,
	}
}

func (lk *LinkToCompressed) String() string {
	return fmt.Sprintf("original sequence id: %d, reference range: (%d, %d)",
		lk.OrgSeqId, lk.RefStart, lk.RefEnd)
}
