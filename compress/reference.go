package main

import (
	"strings"

	"github.com/BurntSushi/cablastp"
)

// referenceDB represents a set of unique sequences that comprise the "coarse"
// database. Sequences in the ReferenceDB are use to re-create the original
// sequences.
type referenceDB struct {
	seqs  []*cablastp.ReferenceSeq
	seeds seeds
}

// newReferenceDB takes a list of initial original sequences, and adds each
// sequence to the reference database unchanged. Seeds are also generated for
// each K-mer in each original sequence.
func newReferenceDB(orgSeqs []*cablastp.OriginalSeq) *referenceDB {
	refdb := &referenceDB{
		seqs:  make([]*cablastp.ReferenceSeq, 0, len(orgSeqs)),
		seeds: newSeeds(),
	}
	for _, orgSeq := range orgSeqs {
		refdb.add(orgSeq)
	}
	return refdb
}

// Add takes an original sequence, converts it to a reference sequence, and
// adds it as a new reference sequence to the reference database. Seeds are
// also generated for each K-mer in the sequence. The resulting reference
// sequence is returned.
func (refdb *referenceDB) add(
	orgSeq *cablastp.OriginalSeq) *cablastp.ReferenceSeq {

	nextIndex := len(refdb.seqs)
	refSeq := cablastp.NewReferenceSeq(nextIndex, orgSeq.Name, orgSeq.Residues)
	refdb.seeds.add(nextIndex, refSeq)
	refdb.seqs = append(refdb.seqs, refSeq)
	return refSeq
}

// String returns a FASTA representation of the reference database.
func (refdb *referenceDB) String() string {
	seqStrs := make([]string, len(refdb.seqs))
	for i, seq := range refdb.seqs {
		seqStrs[i] = seq.String()
	}
	return strings.Join(seqStrs, "\n")
}
