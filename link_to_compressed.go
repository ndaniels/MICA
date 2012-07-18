package cablastp

// LinkToCompressed represents a link from a reference sequence to a
// compressed original sequence. It serves as a bridge from a BLAST hit in
// the coarse database to the corresponding original sequence that is
// redundant to the specified residue range in the reference sequence.
type LinkToCompressed struct {
	OrgSeqId         int
	RefStart, RefEnd int
}

func NewLinkToCompressed(orgSeqId int, refSeq *ReferenceSeq) *LinkToCompressed {
	return &LinkToCompressed{
		OrgSeqId: orgSeqId,
		RefStart: refSeq.Offset,
		RefEnd:   refSeq.Offset + refSeq.Len(),
	}
}
