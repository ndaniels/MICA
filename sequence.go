package main

import (
	"io"

	"github.com/kortschak/biogo/seq"
	"github.com/kortschak/biogo/io/seqio/fasta"
)

type sequence struct {
	*seq.Seq
}

func newSequence(s *seq.Seq) sequence {
	return sequence{
		Seq: s,
	}
}

func readSeqs(fileName string) ([]sequence, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}

	sequences := make([]sequence, 0)
	for {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		sequences = append(sequences, newSequence(seq))
	}

	return sequences, nil
}

