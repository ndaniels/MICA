package main

import (
	"github.com/ndaniels/neutronium"
	"sort"
)

type starterSeq struct {
	oSeq   *neutronium.OriginalSeq
	oSeqId int
}

type BySeqLength []starterSeq

func (s BySeqLength) Len() int {
	return len(s)
}

func (s BySeqLength) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

func (s BySeqLength) Less(i, j int) bool {
	return len(s[i].oSeq.Residues) < len(s[j].oSeq.Residues)
}

func primeCoarseDB(clusterThresh float64, db *neutronium.DB, starterSeqs []starterSeq) {
	skipTable := make(map[int]bool)
	coarsedb := db.CoarseDB
	mem := newMemory()
	sort.Sort(BySeqLength(starterSeqs))
	for rowInd, rowSeq := range starterSeqs {
		if !skipTable[rowInd] {

			comSeq := neutronium.NewCompressedSeq(rowSeq.oSeqId, rowSeq.oSeq.Name)
			addWithoutMatch(&comSeq, coarsedb, rowSeq.oSeqId, rowSeq.oSeq)
			corSeqId := rowInd
			corSeq := coarsedb.CoarseSeqGet(uint(corSeqId))
			corLen := uint(len(corSeq.Residues))

			for colInd, colSeq := range starterSeqs[rowInd:] {
				if !skipTable[colInd] {

					comp := compareSeqs(clusterThresh, corSeqId, colSeq.oSeqId, corSeq, colSeq.oSeq, mem)
					if comp.distance <= clusterThresh {

						skipTable[colInd] = true
						comSeq.Add(neutronium.NewLinkToCoarse(
							uint(corSeqId), 0, corLen, comp.alignment))
						corSeq.AddLink(neutronium.NewLinkToCompressed(
							uint32(colSeq.oSeqId), 0, uint16(corLen)))
					}
				}
			}
		}
	}

}
