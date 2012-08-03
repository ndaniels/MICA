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
func ReadOriginalSeqs(fileName string,
	ignore []byte) (chan ReadOriginalSeq, error) {

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
			for i, residue := range seq.Seq {
				for _, toignore := range ignore {
					if toignore == residue {
						seq.Seq[i] = 'X'
						break
					}
				}
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
func ReadReferenceSeqs(fileName string, f func(seq *CoarseSeq)) error {
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
		f(NewBiogoCoarseSeq(i, seq))
	}
	return nil
}
