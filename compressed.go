package main

import (
	"fmt"
	"io"
	"strings"

	"github.com/kortschak/biogo/io/seqio/fasta"
)

// compressedDb represents the compressed database of reference sequences and
// links from reference sequences to original sequences from the input FASTA
// file.
type compressedDb struct {
	// seqs is a simple list of reference sequences that makes up the core of
	// the compressed database. In essence, this is the data written as output
	// in the "coarse" FASTA file. Note that the indices of this slice are
	// significant, and correspond precisely to indices in the linkTable.
	seqs []*referenceSeq

	// links is a link table, represented by a map indexed by reference
	// sequence indices from 'seqs' to slice of linkEntry.
	links linkTable

	// seeds is a slice of slice of seed locations. The top-level slice is
	// indexed by hashed kmers (where the length of the slice is always
	// alphabetsize ^ kmersize). Each seed location contains a reference
	// sequence index and a reference residue index, so that a seed location,
	// in combination with 'seqs', can be used to pinpoint precisely where
	// the kmer starts in the reference sequence.
	seeds seeds
}

// newCompressedDb creates a new reference database and initializes it to a set
// of initial sequences. It also initializes the seeds table.
func newCompressedDb(seqs []*originalSeq) *compressedDb {
	refSeqs := make([]*referenceSeq, len(seqs))
	for i, seq := range seqs {
		refSeqs[i] = newReferenceSeq(seq.Seq)
	}
	return &compressedDb{
		seqs: refSeqs,
		links: newLinkTable(),
		seeds: newSeeds(refSeqs),
	}
}

// nextRefIndex returns the next index that should be used when adding a
// reference sequence to the compressed database. This value increases by 1
// after a new sequence is added.
func (cdb *compressedDb) nextRefIndex() int {
	return len(cdb.seqs)
}

// add adds an originalSeq sequence to the reference database. This process is
// the meat and potatoes of cablast compression.
//
// An original sequence may result in a *combination of the following things:
// 1) Multiple additions to the compressed database (multiple reference
// sequences).
// 2) The seeds table updated with any additions to the compressed database.
// 3) Multiply linkEntry's added to the linkTable.
func (cdb *compressedDb) add(origSeq *originalSeq) {
	// Keep track of two pointers. 'current' refers to the residue index in the
	// original sequence that extension is currently originating from.
	// 'lastMatch' refers to the residue index of the *end* of the last match
	// with a reference sequence in the compressed database.
	lastMatch, current := 0, 0

	// Iterate through the original sequence a 'kmer' at a time.
	for i := 0; i < len(origSeq.residues()) - flagSeedSize; i++ {
		kmer := origSeq.residues()[i:i+flagSeedSize]
		if !allUpperAlpha(kmer) {
			continue
		}

		seeds := cdb.seeds.lookup(kmer)
		possibleMatches := make([]linkEntry, 0, len(seeds))

		// Screw ungapped extension. Just expand in both directions and run
		// Smith-Waterman. If there's a decent value returned, use that.

		// for _, seedLoc := range seeds { 
			// ??? 
			// possibleMatches = append(possibleMatches, linkEntry) 
		// } 
		if len(possibleMatches) > 0 {
			bestMatch := possibleMatches[0]
			for _, possibleMatch := range possibleMatches[1:] {
				if bestMatch.Less(possibleMatch) {
					bestMatch = possibleMatch
				}
			}

			cdb.addToCompressed(origSeq.subSequence(
				lastMatch, bestMatch.original.origStartRes))

			cdb.links.add(bestMatch)
			current = bestMatch.original.origEndRes
			lastMatch = current
		}

		current++
	}
	cdb.addToCompressed(origSeq.subSequence(lastMatch, origSeq.Len()))
}

func (cdb *compressedDb) addToCompressed(subOrigSeq *originalSeq) {
	refSeq := newReferenceSeq(subOrigSeq.sequence.Seq)
	cdb.seqs = append(cdb.seqs, refSeq)
	cdb.seeds.add(len(cdb.seqs), refSeq)
}

// residues returns the residues in the compressed database for the sequence
// pointed to by seedLoc with offsets start and end, where the offsets apply
// to the residue index pointed to by seedLoc.
//
// Specifically, the residues returned for some 'seq' sequence:
// seq[residueIndex - backOff : residueIndex + fwdOff]
//
// TODO: Bounds checking.
func (cdb *compressedDb) residues(loc seedLoc, backOff, fwdOff int) []byte {
	sequence := cdb.seqs[loc.seqInd].residues()
	start, end := loc.resInd - backOff, loc.resInd + fwdOff
	return sequence[max(0, start):min(len(sequence) - 1, end)]
}

func (cdb *compressedDb) String() string {
	strs := make([]string, len(cdb.seqs))
	for i, seq := range cdb.seqs {
		strs[i] = fmt.Sprintf("> %s\n%s", seq.ID, string(seq.residues()))
	}
	return strings.Join(strs, "\n")
}

// seedLoc represents the information required to translate a seed to a slice
// of residues from the reference database. Namely, the index of the sequence
// in the reference database and the index of the residue where the seed starts
// in that sequence. The length of the seed is a constant set at
// run-time: flagSeedSize.
type seedLoc struct {
	seqInd, resInd int
}

func newSeedLoc(seqInd, resInd int) seedLoc {
	return seedLoc{
		seqInd: seqInd,
		resInd: resInd,
	}
}

func (loc seedLoc) String() string {
	return fmt.Sprintf("(%d, %d)", loc.seqInd, loc.resInd)
}

type seeds [][]seedLoc

func newSeeds(refSeqs []*referenceSeq) seeds {
	seeds := make(seeds, pow(alphaSize, flagSeedSize))
	for seqInd, seq := range refSeqs {
		seeds.add(seqInd, seq)
	}
	return seeds
}

func (ss seeds) add(refSeqIndex int, refSeq *referenceSeq) {
	for i := 0; i < len(refSeq.residues()) - flagSeedSize; i++ {
		kmer := refSeq.residues()[i:i+flagSeedSize]
		if !allUpperAlpha(kmer) {
			continue
		}

		kmerIndex := hashKmer(kmer)
		loc := newSeedLoc(refSeqIndex, i)

		// If no memory has been allocated for this kmer, then do so now.
		if ss[kmerIndex] == nil {
			ss[kmerIndex] = make([]seedLoc, 1)
			ss[kmerIndex][0] = loc
		} else {
			ss[kmerIndex] = append(ss[kmerIndex], loc)
		}
	}
}

func (ss seeds) lookup(kmer []byte) []seedLoc {
	return ss[hashKmer(kmer)]
}

func (ss seeds) String() string {
	strs := make([]string, 0, len(ss))
	for key, locs := range ss {
		if locs == nil {
			continue
		}

		lstrs := make([]string, len(locs))
		for i, loc := range locs {
			lstrs[i] = loc.String()
		}

		strs = append(strs,
			fmt.Sprintf("%d: %s", key, strings.Join(lstrs, " ")))
	}
	return strings.Join(strs, "\n")
}

func readSeqs(fileName string) ([]*originalSeq, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}

	sequences := make([]*originalSeq, 0)
	for i := 0; true; i++ {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		sequences = append(sequences, newOriginalSeq(i, seq))
	}

	return sequences, nil
}
