package main

import(
	"github.com/BurntSushi/cablastp"
)

type starterSeq struct {
	oSeq *cablastp.OriginalSeq
	oSeqId int
}



func initCoarseDB(clusterThresh float64, db *cablastp.DB, starterSeqs []starterSeq){
	skipTable := make( map[int]bool )
	coarsedb := db.CoarseDB
	mem := newMemory()

	for rowInd, rowSeq := starterSeqs {
		if ! skipTable[rowInd]{

			comSeq := cablastp.NewCompressedSeq(rowSeq.oSeqId, rowSeq.oSeq.Name)
			addWithoutMatch(&comSeq, coarsedb, rowSeq.oSeqId, rowSeq.oSeq)
			corSeqId := rowInd
			corSeq := coarsedb.CoarseSeqGet( uint(corSeqId))
			corLen := uint( len( corSeq.Residues))

			for colInd, colSeq := starterSeqs[rowInd:] {
				if ! skipTable[colInd]{

					comp := compareSeqs(clusterThresh, comSeq, colSeq.oSeq, mem)
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