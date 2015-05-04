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

func primeCoarseDB(maxRadius float64, db *neutronium.DB, seedTable *neutronium.SeedTable, starterSeqs []starterSeq) {
	primeProgressBar := neutronium.ProgressBar{
		Label:   "Priming Database",
		Total:   uint64(len(starterSeqs)),
		Current: uint64(0),
	}

	skipTable := make(map[int]bool)
	coarsedb := db.CoarseDB
	mem := newMemory()
	sort.Sort(BySeqLength(starterSeqs))
	for rowInd, rowSeq := range starterSeqs {
		primeProgressBar.ClearAndDisplay()
		if !skipTable[rowInd] {
			comRowSeq := neutronium.NewCompressedSeq(rowSeq.oSeqId, rowSeq.oSeq.Name)
			corSeqId := addWithoutMatch(&comRowSeq, coarsedb, rowSeq.oSeqId, rowSeq.oSeq, seedTable)
			corSeq := coarsedb.CoarseSeqGet(uint(corSeqId))
			// corLen := uint(len(corSeq.Residues))

			for colPos, colSeq := range starterSeqs[rowInd:] {
				colInd := colPos + rowInd + 1
				if !skipTable[colInd] {

					comp := compareSeqs(maxRadius, corSeqId, colSeq.oSeqId, corSeq, colSeq.oSeq, seedTable, mem)
					//neutronium.Vprintf("%d distance of %f \n", colInd, comp.distance)
					if comp.distance <= maxRadius {

						skipTable[colInd] = true

						addWithMatch(colSeq.oSeq, corSeq, comp.alignment, colSeq.oSeqId, corSeqId)

						// comColSeq := neutronium.NewCompressedSeq(colSeq.oSeqId, colSeq.oSeq.Name)
						// comColSeq.Add(neutronium.NewLinkToCoarse(
						// 	uint(corSeqId), 0, corLen, comp.alignment))
						// corSeq.AddLink(neutronium.NewLinkToCompressed(
						// 	uint32(colSeq.oSeqId), 0, uint16(corLen)))
					}
				}
			}
		}
		primeProgressBar.Increment()
	}

}
