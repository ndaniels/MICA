package cablastp

type CompressedDB []*CompressedSeq

func (comdb *CompressedDB) Add(comSeq *CompressedSeq) {
	comdb = append(comdb, comSeq)
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


