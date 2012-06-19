package main

import (
	"github.com/kortschak/biogo/seq"
)

// linkTable corresponds to a table of link entries, where rows correspond
// to reference sequences, and columns correspond to link entries (which each
// contain a pointer to an original location).
type linkTable map[int][]linkEntry

func newLinkTable() linkTable {
	return make(linkTable, 100)
}

func (lt linkTable) add(lkEntry linkEntry) {
	ind := lkEntry.refSeqId
	if _, ok := lt[ind]; !ok {
		lt[ind] = make([]linkEntry, 1)
		lt[ind][0] = lkEntry
	} else {
		lt[ind] = append(lt[ind], lkEntry)
	}
}

type linkEntry struct {
	refSeqId int
	refStartRes, refEndRes int
	original originalLoc
}

func newLinkEntry(refSeqId, refStart, refEnd int, origSeq *originalSeq,
	align seq.Alignment) linkEntry {

	return linkEntry{
		refSeqId: refSeqId,
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

