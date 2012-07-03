package main

import (
	"github.com/kortschak/biogo/seq"
)

// linkTable corresponds to a table of link entries, where rows correspond
// to reference sequences, and columns correspond to link entries (which each
// contain a pointer to an original location).
type linkTable map[int][]linkEntry

// newLinkTabe constructs a new map for a link table with some initial memory.
// The linkEntry slices themselves are not allocated, though. (That's done in
// 'add' when necessary.)
func newLinkTable() linkTable {
	return make(linkTable, 100)
}

// add adds a linkEntry to the given linkTable, with the index of the reference
// sequence in the compressed database acting as the index in the link table.
// If this is the first link entry for that reference sequence, a new slice
// of linkEntry is allocated. Otherwise, the linkEntry is appended to the end
// of the existing slice.
func (lt linkTable) add(lkEntry linkEntry) {
	ind := lkEntry.refSeqId
	if _, ok := lt[ind]; !ok {
		lt[ind] = make([]linkEntry, 1)
		lt[ind][0] = lkEntry
	} else {
		lt[ind] = append(lt[ind], lkEntry)
	}
}

// linkEntry represents a link between a part of a reference sequence already 
// in the compressed database, and a part of an original sequence from the
// input FASTA file.
//
// Links are created whenever a match (for some definition of 'match') between
// a part of an original sequence and a part of a reference sequence is made.
type linkEntry struct {
	// refSeqId is the index of the reference sequence in the compressed
	// database. It will also be the index of the row of linkEntry in the
	// linkTable.
	refSeqId int

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
func newLinkEntry(refSeqId, refStart, refEnd int, origSeq *originalSeq,
	align seq.Alignment) linkEntry {

	return linkEntry{
		refSeqId: refSeqId,
		refStartRes: refStart,
		refEndRes: refEnd,
		original: newOriginalLoc(origSeq, align),
	}
}

// Less implements an ordinal relationship between two link entries.
// Currently, one linkEntry is "better" than the other if it corresponds to a
// bigger part of a reference sequence in the compressed database.
func (le1 linkEntry) Less(le2 linkEntry) bool {
	return (le1.refEndRes - le1.refStartRes) < (le2.refEndRes - le2.refStartRes)
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
		diff: newEditScript(align),
		origSeqId: origSeq.seqId,
		origStartRes: origSeq.start,
		origEndRes: origSeq.end,
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

