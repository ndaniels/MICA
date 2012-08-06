package cablastp

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"

	"code.google.com/p/biogo/io/seqio/fasta"
)

func (coarsedb *CoarseDB) readFasta() error {
	fastaReader := fasta.NewReader(coarsedb.FileFasta)
	for i := 0; true; i++ {
		seq, err := fastaReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		coarsedb.Seqs = append(coarsedb.Seqs, NewBiogoCoarseSeq(i, seq))
	}
	return nil
}

func (coarsedb *CoarseDB) saveFasta() error {
	bufWriter := bufio.NewWriter(coarsedb.FileFasta)
	for i := coarsedb.seqsRead; i < len(coarsedb.Seqs); i++ {
		_, err := fmt.Fprintf(bufWriter,
			"> %d\n%s\n", i, string(coarsedb.Seqs[i].Residues))
		if err != nil {
			return err
		}
	}
	if err := bufWriter.Flush(); err != nil {
		return err
	}
	return nil
}

func (coarsedb *CoarseDB) readSeeds() error {
	var err, err2 error

	gzipReader, err := gzip.NewReader(coarsedb.FileSeeds)
	if err != nil {
		return err
	}
	bufReader := bufio.NewReader(gzipReader)

	var byteReader *bytes.Reader
	var hash, seqInd int
	var resInd int16
	for {
		line, err := bufReader.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		byteReader = bytes.NewReader(line[:len(line)-1])
		err = binary.Read(byteReader, binary.BigEndian, &hash)
		if err != nil {
			return err
		}
		for {
			err = binary.Read(byteReader, binary.BigEndian, &seqInd)
			err2 = binary.Read(byteReader, binary.BigEndian, &resInd)
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err2 != nil {
				return err2
			}

			sl := NewSeedLoc(seqInd, resInd)
			if coarsedb.Seeds.Locs[hash] == nil {
				coarsedb.Seeds.Locs[hash] = sl
			} else {
				lk := coarsedb.Seeds.Locs[hash]
				for ; lk.Next != nil; lk = lk.Next {
				}
				lk.Next = sl
			}
		}
	}
	if err := gzipReader.Close(); err != nil {
		return err
	}
	return nil
}

func (coarsedb *CoarseDB) saveSeeds() error {
	var i int32

	gzipWriter, err := gzip.NewWriterLevel(coarsedb.FileSeeds, gzip.BestSpeed)
	if err != nil {
		return err
	}
	for i = 0; i < int32(coarsedb.Seeds.powers[coarsedb.Seeds.SeedSize]); i++ {
		if coarsedb.Seeds.Locs[i] == nil {
			continue
		}

		binary.Write(gzipWriter, binary.BigEndian, i)
		for loc := coarsedb.Seeds.Locs[i]; loc != nil; loc = loc.Next {
			binary.Write(gzipWriter, binary.BigEndian, loc.SeqInd)
			binary.Write(gzipWriter, binary.BigEndian, loc.ResInd)
		}
		if _, err := gzipWriter.Write([]byte{'\n'}); err != nil {
			return err
		}
	}
	if err := gzipWriter.Close(); err != nil {
		return err
	}
	return nil
}

func (coarsedb *CoarseDB) saveSeedsPlain() error {
	csvWriter := csv.NewWriter(coarsedb.plainSeeds)
	record := make([]string, 0, 10)
	for i := 0; i < coarsedb.Seeds.powers[coarsedb.Seeds.SeedSize]; i++ {
		if coarsedb.Seeds.Locs[i] == nil {
			continue
		}

		record = record[:0]
		record = append(record, string(coarsedb.Seeds.unhashKmer(i)))
		for loc := coarsedb.Seeds.Locs[i]; loc != nil; loc = loc.Next {
			record = append(record,
				fmt.Sprintf("%d", loc.SeqInd),
				fmt.Sprintf("%d", loc.ResInd))
		}
		if err := csvWriter.Write(record); err != nil {
			return err
		}
	}
	csvWriter.Flush()
	return nil
}

func (coarsedb *CoarseDB) readLinks() error {
	var err, err2, err3 error
	gzipReader, err := gzip.NewReader(coarsedb.FileLinks)
	if err != nil {
		return err
	}
	bufReader := bufio.NewReader(gzipReader)

	var byteReader *bytes.Reader
	var orgSeqId int
	var coarseStart, coarseEnd int16

	for coarseSeqId := 0; true; coarseSeqId++ {
		line, err := bufReader.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		byteReader = bytes.NewReader(line[:len(line)-1])
		for {
			err = binary.Read(byteReader, binary.BigEndian, &orgSeqId)
			err2 = binary.Read(byteReader, binary.BigEndian, &coarseStart)
			err3 = binary.Read(byteReader, binary.BigEndian, &coarseEnd)
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err2 != nil {
				return err2
			}
			if err3 != nil {
				return err3
			}

			newLink := NewLinkToCompressed(orgSeqId, coarseStart, coarseEnd)
			coarsedb.Seqs[coarseSeqId].addLink(newLink)
		}
	}
	if err := gzipReader.Close(); err != nil {
		return err
	}
	return nil
}

func (coarsedb *CoarseDB) saveLinks() error {
	gzipWriter, err := gzip.NewWriterLevel(coarsedb.FileLinks, gzip.BestSpeed)
	if err != nil {
		return err
	}
	for _, seq := range coarsedb.Seqs {
		for link := seq.Links; link != nil; link = link.Next {
			binary.Write(gzipWriter, binary.BigEndian, link.OrgSeqId)
			binary.Write(gzipWriter, binary.BigEndian, link.CoarseStart)
			binary.Write(gzipWriter, binary.BigEndian, link.CoarseEnd)
		}
		if _, err := gzipWriter.Write([]byte{'\n'}); err != nil {
			return err
		}
	}
	if err := gzipWriter.Close(); err != nil {
		return nil
	}
	return nil
}

func (coarsedb *CoarseDB) saveLinksPlain() error {
	csvWriter := csv.NewWriter(coarsedb.plainLinks)
	record := make([]string, 0, 10)
	for _, seq := range coarsedb.Seqs {
		record = record[:0]
		for link := seq.Links; link != nil; link = link.Next {
			record = append(record,
				fmt.Sprintf("%d", link.OrgSeqId),
				fmt.Sprintf("%d", link.CoarseStart),
				fmt.Sprintf("%d", link.CoarseEnd))
		}
		if err := csvWriter.Write(record); err != nil {
			return err
		}
	}
	csvWriter.Flush()
	return nil
}

func (comdb *CompressedDB) ReadSeq(
	coarsedb *CoarseDB, orgSeqId int) (OriginalSeq, error) {

	off, err := comdb.orgSeqOffset(orgSeqId)
	if err != nil {
		return OriginalSeq{}, err
	}

	newOff, err := comdb.File.Seek(off, 0)
	if err != nil {
		return OriginalSeq{}, err
	} else if newOff != off {
		return OriginalSeq{},
			fmt.Errorf("Tried to seek to offset %d in the compressed "+
				"database, but seeked to %d instead.", off, newOff)
	}

	record, err := comdb.csvReader.Read()
	if err != nil {
		return OriginalSeq{}, err
	}

	cseq, err := readCompressedSeq(orgSeqId, record)
	if err != nil {
		return OriginalSeq{}, err
	}
	return cseq.Decompress(coarsedb)
}

func (comdb *CompressedDB) orgSeqOffset(id int) (seqOff int64, err error) {
	tryOff := int64(id) * 8
	realOff, err := comdb.Index.Seek(tryOff, 0)
	if err != nil {
		return 0, err
	} else if tryOff != realOff {
		return 0,
			fmt.Errorf("Tried to seek to offset %d in the compressed index, "+
				"but seeked to %d instead.", tryOff, realOff)
	}
	err = binary.Read(comdb.Index, binary.BigEndian, &seqOff)
	return
}

func readCompressedSeq(id int, record []string) (CompressedSeq, error) {
	cseq := CompressedSeq{
		Id:    id,
		Name:  string([]byte(record[0])),
		Links: make([]LinkToCoarse, (len(record)-1)/4),
	}

	for i := 1; i < len(record); i += 4 {
		coarseSeqId64, err := strconv.Atoi(record[i+0])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseStart64, err := strconv.Atoi(record[i+1])
		if err != nil {
			return CompressedSeq{}, nil
		}
		coarseEnd64, err := strconv.Atoi(record[i+2])
		if err != nil {
			return CompressedSeq{}, nil
		}
		lk := NewLinkToCoarseNoDiff(
			int(coarseSeqId64), int(coarseStart64), int(coarseEnd64))
		lk.Diff = string([]byte(record[i+3]))

		cseq.Add(lk)
	}
	return cseq, nil
}

func (comdb *CompressedDB) writer() {
	var record []string
	var err error

	byteOffset := int64(0)
	buf := new(bytes.Buffer)
	csvWriter := csv.NewWriter(buf)
	csvWriter.Comma = ','
	csvWriter.UseCRLF = false

	for cseq := range comdb.writerChan {
		// Reset the buffer so it's empty. We want it to only contain
		// the next record we're writing.
		buf.Reset()

		// Allocate memory for creating the next record.
		// A record is a sequence name followed by four-tuples of links:
		// (coarse-seq-id, coarse-start, coarse-end, diff).
		record = make([]string, 0, 1+4*len(cseq.Links))
		record = append(record, cseq.Name)
		for _, link := range cseq.Links {
			record = append(record,
				fmt.Sprintf("%d", link.CoarseSeqId),
				fmt.Sprintf("%d", link.CoarseStart),
				fmt.Sprintf("%d", link.CoarseEnd),
				link.Diff)
		}

		// Write the record to our *buffer* and flush it.
		if err = csvWriter.Write(record); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}
		csvWriter.Flush()

		// Pass the bytes on to the compressed file.
		if _, err = comdb.File.Write(buf.Bytes()); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}

		// Now write the byte offset that points to the start of this record.
		err = binary.Write(comdb.Index, binary.BigEndian, byteOffset)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			os.Exit(1)
		}

		// Increment the byte offset to be at the end of this record.
		byteOffset += int64(buf.Len())
	}
	comdb.File.Close()
	comdb.writerDone <- struct{}{}
}
