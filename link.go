package main

import (
	"github.com/kortschak/biogo/seq"
)

// linkTable corresponds to a table of link entries, where rows correspond
// to reference sequences, and columns correspond to link entries (which each
// contain a pointer to an original location).
type linkTable [][]linkEntry

func newLinkTable() *linkTable {
	lt := make(linkTable, 0, 100)
	return &lt
}

func (lt *linkTable) add(refSeqIndex, lkEntry linkEntry) {
	if len(lt[refSeqIndex]) == 0 {
		lt[refSeqIndex] = make([]linkEntry, 1)
		lt[refSeqIndex][0] = lkEntry
	} else {
		lt[refSeqIndex] = append(lt[refSeqIndex], lkEntry)
	}
}

type linkEntry struct {
	refStartRes, refEndRes int
	original originalLoc
}

func newLinkEntry(refStart, refEnd int, origSeq *originalSeq,
	align seq.Alignment) linkEntry {

	return linkEntry{
		refStartRes: refStart,
		refEndRes: refEnd,
		original: newOriginalLoc(origSeq, align),
	}
}

func (le1 linkEntry) Less(le2 linkEntry) bool {
	return (le1.refEndRes - le1.refStartRes) < (le2.refEndRes - le2.refStartRes)
}

type originalLoc struct {
	diff editScript
	origSeqId int
	origStartRes int
	origEndRes int
}

func newOriginalLoc(origSeq *originalSeq, align seq.Alignment) originalLoc {
	return originalLoc{
		diff: newEditScript(align),
		origSeqId: origSeq.seqId,
		origStartRes: origSeq.start,
		origEndRes: origSeq.end,
	}
}

type editScript string

func newEditScript(align seq.Alignment) editScript {
	return ""
}

