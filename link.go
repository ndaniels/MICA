package main

import (
	"github.com/kortschak/biogo/seq"
)

type linkTable [][]linkEntry

func newLinkTable() linkTable {
	return make(linkTable, 0, 100)
}

func (lt linkTable) add() {
}

type linkEntry struct {
	refStartRes, refEndRes int
	original originalLoc
}

type originalLoc struct {
	diff editScript
	origSeqId int
	origStartRes int
}

type editScript string

func newEditScript(align seq.Alignment) {
}

