package main

import (
	"github.com/kortschak/biogo/seq"
)

// linkEntry represents a link between a part of a reference sequence already 
// in the compressed database, and a part of an original sequence from the
// input FASTA file.
//
// Links are created whenever a match (for some definition of 'match') between
// a part of an original sequence and a part of a reference sequence is made.
type linkEntry struct {
	// refStartRes and refEndRes correspond to amino acid residues indexes
	// into the reference sequence's 'residues' member.
	refStartRes, refEndRes int

	// original encapsulates all information about the original input sequence,
	// including the start and end residues, the edit script, and the index
	// of the original sequence in the input FASTA file.
	original originalLoc
}

// newLinkEntry creates a new linkEntry value suitable to add to the linkTable.
// Note that an originalLoc value is created for you with the given information.
func newLinkEntry(refStart, refEnd int, origSeq *originalSeq,
	align seq.Alignment) *linkEntry {

	return &linkEntry{
		refStartRes: refStart,
		refEndRes:   refEnd,
		original:    newOriginalLoc(origSeq, align),
	}
}

// originalLoc encapsulates the information related to an original sequence in
// each linkEntry. Of particular note is the edit script, which when "applied"
// to the corresponding reference sequence in the linkEntry, will return
// precisely the amino acid residue sequence of the original sequence.
type originalLoc struct {
	// diff represents the difference between this original sequence and the
	// reference sequence. It serves as a compressed mechanism of retrieving
	// the original sequence from the reference sequence.
	diff editScript

	// origSeqId is the index of the original sequence in the input FASTA file.
	// It is analogous to refSeqId in linkEntry.
	origSeqId int

	// origStartRes and origEndRes represent the window of the original sequence
	// that is mapped to a particular piece of a reference sequence in the
	// compressed database.
	origStartRes, origEndRes int
}

// newOriginalLoc creates a new originalLoc value from an original sequence
// and an alignment with a piece of a reference sequence from the compressed
// database.
func newOriginalLoc(origSeq *originalSeq, align seq.Alignment) originalLoc {
	return originalLoc{
		diff:         newEditScript(align),
		origSeqId:    origSeq.id,
		origStartRes: origSeq.offset,
		origEndRes:   origSeq.offset + origSeq.Len(),
	}
}

// editScript corresponds to a small 'diff' *from* a reference sequence *to*
// an original sequence. When applied to a reference sequence, the original
// sequence is returned.
type editScript string

// newEditScript writes an edit script from an alignment.
func newEditScript(align seq.Alignment) editScript {
	return ""
}
