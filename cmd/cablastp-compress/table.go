package main

import(
	"github.com/BurntSushi/cablastp"
	"sort"
)

type starterSeq struct {
	oSeq *cablastp.OriginalSeq
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



func primeCoarseDB(clusterThresh float64, db *cablastp.DB, starterSeqs []starterSeq){
	skipTable := make( map[int]bool )
	coarsedb := db.CoarseDB
	mem := newMemory()

	sort.Sort( BySeqLength(starterSeqs))

	for rowInd, rowSeq := starterSeqs {
		if ! skipTable[rowInd]{

			comSeq := cablastp.NewCompressedSeq(rowSeq.oSeqId, rowSeq.oSeq.Name)
			addWithoutMatch(&comSeq, coarsedb, rowSeq.oSeqId, rowSeq.oSeq)
			corSeqId := rowInd
			corSeq := coarsedb.CoarseSeqGet( uint(corSeqId))
			corLen := uint( len( corSeq.Residues))

			for colInd, colSeq := starterSeqs[rowInd:] {
				if ! skipTable[colInd]{

					comp := compareSeqs(clusterThresh, corSeq, colSeq.oSeq, mem)
					if comp.distance <= clusterThresh {

						skipTable[colInd] = true
						comSeq.Add(cablastp.NewLinkToCoarse(
									uint(corSeqId), 0, corLen, comp.alignment))
						corSeq.AddLink(cablastp.NewLinkToCompressed(
									uint32(colSeq.oSeqId), 0, corLen))
					}
				}
			}
		}
	}

}