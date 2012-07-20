package cablastp

import (
	"io"

	"code.google.com/p/biogo/io/seqio/fasta"
)

type ReadOriginalSeq struct {
	Seq *OriginalSeq
	Err error
}

// ReadOriginalSeqs reads a FASTA formatted file and returns a slice of
// original sequences.
func ReadOriginalSeqs(fileName string) (chan ReadOriginalSeq, error) {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return nil, err
	}
	seqChan := make(chan ReadOriginalSeq, 200)
	go func() {
		for i := 0; true; i++ {
			seq, err := reader.Read()
			if err == io.EOF {
				close(seqChan)
				break
			}
			if err != nil {
				seqChan <- ReadOriginalSeq{
					Seq: nil,
					Err: err,
				}
				close(seqChan)
				break
			}
			seqChan <- ReadOriginalSeq{
				Seq: NewBiogoOriginalSeq(i, seq),
				Err: nil,
			}
		}
	}()
	return seqChan, nil
}

// ReadReferenceSeqs reads a FASTA formatted file and returns a slice of
// reference sequences.
func ReadReferenceSeqs(fileName string, f func(seq *ReferenceSeq)) error {
	reader, err := fasta.NewReaderName(fileName)
	if err != nil {
		return err
	}
	for i := 0; true; i++ {
		seq, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		f(NewBiogoReferenceSeq(i, seq))
	}
	return nil
}
