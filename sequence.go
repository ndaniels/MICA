package main

import (
	"fmt"
	"io"
	"strings"

	"github.com/kortschak/biogo/io/seqio/fasta"
	"github.com/kortschak/biogo/seq"
)

type reference []referenceSeq

// newReference creates a new reference database and initializes it to a set
// of initial sequences. It also initializes the seeds table.
func newReference(seqs []originalSeq) (reference, seeds) {
	ref := make(reference, len(seqs))
	for i, seq := range seqs {
		ref[i] = newReferenceSeq(seq.Seq)
	}

	seeds := newSeeds(ref)
	seeds.add(ref)
	return ref, seeds
}

func (refs reference) String() string {
	strs := make([]string, len(refs))
	for i, seq := range refs {
		strs[i] = fmt.Sprintf("> %s\n%s", seq.ID, string(seq.residues()))
	}
	return strings.Join(strs, "\n")
}

type referenceSeq struct {
	*seq.Seq
}

func newReferenceSeq(s *seq.Seq) referenceSeq {
	s.Seq = []byte(strings.ToUpper(string(s.Seq)))
	return referenceSeq{s}
}

func (seq referenceSeq) residues() []byte {
	return seq.Seq.Seq
}

type location struct {
	seqInd, resInd int
}

func newLocation(seqInd, resInd int) location {
	return location{
		seqInd: seqInd,
		resInd: resInd,
	}
}

func (loc location) String() string {
	return fmt.Sprintf("(%d, %d)", loc.seqInd, loc.resInd)
}

type seeds [][]location

func newSeeds(ref reference) seeds {
	return make(seeds, pow(alphaSize, flagSeedSize))
}

func (ss seeds) add(refs []referenceSeq) {
	for seqInd, seq := range refs {
		for i := 0; i < len(seq.residues()) - flagSeedSize; i++ {
			kmer := seq.residues()[i:i+flagSeedSize]
			kmerIndex := hashKmer(kmer)
			loc := newLocation(seqInd, i)

			// If no memory has been allocated for this kmer, then do so now.
			if ss[kmerIndex] == nil {
				ss[kmerIndex] = make([]location, 1)
				ss[kmerIndex][0] = loc
			} else {
				ss[kmerIndex] = append(ss[kmerIndex], loc)
			}
		}
	}
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

type originalSeq struct {
	*seq.Seq
}

func newOriginalSeq(s *seq.Seq) originalSeq {
	s.Seq = []byte(strings.ToUpper(string(s.Seq)))
	return originalSeq{s}
}

func readSeqs(fileName string) ([]originalSeq, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}

	sequences := make([]originalSeq, 0)
	for {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		sequences = append(sequences, newOriginalSeq(seq))
	}

	return sequences, nil
}
