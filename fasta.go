package cablastp

import (
	"io"

	"code.google.com/p/biogo/io/seqio/fasta"
)

// ReadOriginalSeq is the value sent over `chan ReadOriginalSeq` when a new
// sequence is read from a fasta file.
type ReadOriginalSeq struct {
	Seq *OriginalSeq
	Err error
}

// ReadOriginalSeqs reads a FASTA formatted file and returns a channel that
// each new sequence is sent to.
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
