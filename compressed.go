package cablastp

type CompressedDB struct {
	Seqs []*CompressedSeq
}

func NewCompressedDB() *CompressedDB {
	return &CompressedDB{
		Seqs: make([]*CompressedSeq, 0, 100),
	}
}

func (comdb *CompressedDB) Add(comSeq *CompressedSeq) {
	comdb.Seqs = append(comdb.Seqs, comSeq)
}

// CompressedSeq corresponds to the components of a compressed sequence.
type CompressedSeq struct {
	// Name is an uncompressed string from the original FASTA header.
	Name string

	// Links is an ordered lists of links to portions of the reference
	// database. When all links are followed, the concatenation of each
	// sequence correspond to each link equals the entire original sequence.
	Links []*LinkToReference
}


