package main

import (
	"fmt"
	"io"
	"strings"

	"code.google.com/p/biogo/io/seqio/fasta"
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
		refSeqs[i] = newReferenceSeq(i, seq.name, seq.residues)
	}
	return &compressedDb{
		seqs:  refSeqs,
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
// An original sequence may result in a combination of the following things:
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
	for current = 0; current < origSeq.Len()-flagSeedSize; current++ {
		kmer := origSeq.residues[current : current+flagSeedSize]
		if !allUpperAlpha(kmer) {
			continue
		}

		seeds := cdb.seeds.lookup(kmer)
		// fmt.Println("SEED INFO", current, string(kmer), len(seeds)) 
		possibleMatches := make([]match, 0, len(seeds))

		for _, seedLoc := range seeds {
			// fmt.Println(current, seedLoc.seqInd, seedLoc.resInd) 
			// we need a reference sequence and an original sequence.
			// They may both be subsequences.

			// 1) Attempt ungapped extension. This returns some offset,
			// which may be zero.
			// 2) Attempt gapped extension from the previous offset up to
			// GAPPED_WINDOW_LENGTH.
			// 3) Repeat until gapped extension returns an alignment with
			// less than the sequence identity threshold.
			rseq := cdb.seqs[seedLoc.seqInd]
			startRseq := rseq.newSubSequence(seedLoc.resInd, rseq.Len())
			startOseq := origSeq.newSubSequence(current, origSeq.Len())
			matchPos := 0
			for {
				if matchPos == startOseq.Len() {
					break
				}

				subRseq := startRseq.newSubSequence(matchPos, startRseq.Len())
				subOseq := startOseq.newSubSequence(matchPos, startOseq.Len())

				matchLen := alignUngapped(subRseq, subOseq)
				matchPos += matchLen

				tmpRseq := startRseq.newSubSequence(
					matchPos,
					min(startRseq.Len(), matchPos+flagGappedWindowSize))
				tmpOseq := startOseq.newSubSequence(
					matchPos,
					min(startOseq.Len(), matchPos+flagGappedWindowSize))
				alignment := alignGapped(tmpRseq, tmpOseq)
				id := identity(alignment[0].Seq, alignment[1].Seq)
				if id < flagSeqIdThreshold {
					break
				}

				matchPos += tmpOseq.Len()
			}

			if matchPos-current >= flagMinMatchLen {
				subRseq := rseq.newSubSequence(seedLoc.resInd, matchPos)
				subOseq := origSeq.newSubSequence(current, matchPos)
				alignment := alignGapped(subRseq, subOseq)

				fmt.Println("current to matchPos", current, matchPos)
				fmt.Println("identity",
					identity(alignment[0].Seq, alignment[1].Seq))
				fmt.Println("> ", rseq.name)
				fmt.Println(string(subRseq.residues))
				fmt.Println("> ", origSeq.name)
				fmt.Println(string(subOseq.residues))
				fmt.Println("")
				fmt.Println(alignment)
				fmt.Println("--------------------------------------------")

				link := newLinkEntry(seedLoc.resInd, matchPos, subOseq,
					alignment)
				possibleMatches = append(possibleMatches,
					match{
						rseq: rseq,
						link: link,
					})
			}
		}
		if len(possibleMatches) > 0 {
			bestMatch := possibleMatches[0]
			for _, possibleMatch := range possibleMatches[1:] {
				if bestMatch.Less(possibleMatch) {
					bestMatch = possibleMatch
				}
			}

			if bestMatch.link.original.origStartRes-lastMatch > 0 {
				sub := origSeq.newSubSequence(
					lastMatch, bestMatch.link.original.origStartRes)
				fmt.Println(strings.Repeat("#", 45))
				fmt.Println(sub)
				fmt.Println(strings.Repeat("#", 45))
				cdb.addToCompressed(sub)
			}

			bestMatch.rseq.addLink(bestMatch.link)
			current = bestMatch.link.original.origEndRes
			lastMatch = current
		}
	}
	sub := origSeq.newSubSequence(lastMatch, origSeq.Len())
	fmt.Println(strings.Repeat("#", 45))
	fmt.Println(sub)
	fmt.Println(strings.Repeat("#", 45))
	cdb.addToCompressed(sub)
}

func (cdb *compressedDb) addToCompressed(subOrigSeq *originalSeq) {
	refSeq := newReferenceSeq(cdb.nextRefIndex(),
		subOrigSeq.name, subOrigSeq.residues)
	cdb.seeds.add(len(cdb.seqs), refSeq)
	cdb.seqs = append(cdb.seqs, refSeq)
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
	sequence := cdb.seqs[loc.seqInd].residues
	start, end := loc.resInd-backOff, loc.resInd+fwdOff
	return sequence[max(0, start):min(len(sequence)-1, end)]
}

func (cdb *compressedDb) String() string {
	strs := make([]string, len(cdb.seqs))
	for i, seq := range cdb.seqs {
		strs[i] = seq.String()
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
	for i := 0; i < refSeq.Len()-flagSeedSize; i++ {
		kmer := refSeq.residues[i : i+flagSeedSize]
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
		sequences = append(sequences, newBiogoOriginalSeq(i, seq))
	}

	return sequences, nil
}
