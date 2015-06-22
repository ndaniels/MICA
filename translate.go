package mica

import (
	"bytes"
	// "compress/gzip"
	"io"
	// "os"
	"fmt"
	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/seq"
	"strings"
)

type SearchOperator func(*bytes.Reader) (*bytes.Reader, error)

func TranslateQuerySeqs(
	query *bytes.Reader, action SearchOperator) (*bytes.Reader, error) {

	buf := new(bytes.Buffer)
	f := fasta.NewWriter(buf)
	reader := fasta.NewReader(query)
	for i := 0; true; i++ {
		sequence, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		origSeq := sequence.Bytes()
		n := sequence.Name
		// generate 6 ORFs
		transSeqs := Translate(origSeq)
		for _, s := range transSeqs {
			result := seq.NewSequenceString(n, string(Reduce(s)))
			f.Write(result)
		}

	}

	return bytes.NewReader(buf.Bytes()), nil
}

func Translate(sequence []byte) [][]byte {
	l := len(sequence)
	results := make([][]byte, 0, 6)
	// three ORFs
	for orf := 0; orf <= 2; orf++ {
		var result []byte
		// forward direction
		for i := orf; i < (l - 2); i += 3 {
			var codon []byte
			codon = sequence[i : i+3]
			trans := translate1(codon)
			if trans == '_' {
				// elide stop codons for now?
				continue
			}
			result = append(result, trans)
		}
		results = append(results, result)
		// reverse complement
		result = make([]byte, 0)
		for i := (l - 3 - orf); i >= 0; i -= 3 {
			codon := make([]byte, 3)
			rCodon := sequence[i : i+3]
			codon[0] = complement(rCodon[2])
			codon[1] = complement(rCodon[1])
			codon[2] = complement(rCodon[0])
			trans := translate1(codon)
			if trans == '_' {
				// elide stop codons
				continue
			}
			result = append(result, trans)
		}
		results = append(results, result)
	}
	return results
}

func translate1(codon []byte) byte {
	// if code contains an 'N' -> 'X'
	if strings.ContainsRune(string(codon), 'N') {
		return 'X'
	}
	// otherwise, look up in hash

	var genCode = map[string]byte{
		"ATA": 'I', "ATC": 'I', "ATT": 'I', "ATG": 'M',
		"ACA": 'T', "ACC": 'T', "ACG": 'T', "ACT": 'T',
		"AAC": 'N', "AAT": 'N', "AAA": 'K', "AAG": 'K',
		"AGC": 'S', "AGT": 'S', "AGA": 'R', "AGG": 'R',
		"CTA": 'L', "CTC": 'L', "CTG": 'L', "CTT": 'L',
		"CCA": 'P', "CCC": 'P', "CCG": 'P', "CCT": 'P',
		"CAC": 'H', "CAT": 'H', "CAA": 'Q', "CAG": 'Q',
		"CGA": 'R', "CGC": 'R', "CGG": 'R', "CGT": 'R',
		"GTA": 'V', "GTC": 'V', "GTG": 'V', "GTT": 'V',
		"GCA": 'A', "GCC": 'A', "GCG": 'A', "GCT": 'A',
		"GAC": 'D', "GAT": 'D', "GAA": 'E', "GAG": 'E',
		"GGA": 'G', "GGC": 'G', "GGG": 'G', "GGT": 'G',
		"TCA": 'S', "TCC": 'S', "TCG": 'S', "TCT": 'S',
		"TTC": 'F', "TTT": 'F', "TTA": 'L', "TTG": 'L',
		"TAC": 'Y', "TAT": 'Y', "TAA": '_', "TAG": '_',
		"TGC": 'C', "TGT": 'C', "TGA": '_', "TGG": 'W',
	}
	return genCode[string(codon)]
}

func complement(char byte) byte {
	switch char {
	case 'A':
		return 'T'
	case 'T':
		return 'A'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'N':
		return 'N'
	}
	panic(fmt.Sprintf("bad letter: %c", char))
}
