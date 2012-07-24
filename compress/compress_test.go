package main

import (
	"strings"
	"testing"

	"github.com/BurntSushi/cablastp"
	"github.com/BurntSushi/cablastp/blosum"

	"code.google.com/p/biogo/align/nw"
	"code.google.com/p/biogo/seq"
)

func TestNeedlemanWunsch(t *testing.T) {
	type test struct {
		seq1, seq2 string
		out1, out2 string
	}

	tests := []test{
		{
			"ABCD",
			"ABCD",
			"ABCD",
			"ABCD",
		},
		{
			"GHIKLMNPQR",
			"GAAAHIKLMN",
			"---GHIKLMNPQR",
			"GAAAHIKLMN---",
		},
		{
			"GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
			"---GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
		},
		// { 
		// "ABCDEFGWXYZ", 
		// "ABCDEFMNPQRSTZABEGWXYZ", 
		// "ABCDEF-----------GWXYZ", 
		// "ABCDEFMNPQRSTZABEGWXYZ", 
		// }, 
	}
	sep := strings.Repeat("-", 45)
	for _, test := range tests {
		alignment := nwAlign([]byte(test.seq1), []byte(test.seq2),
			make([]byte, 0, 1000), make([]byte, 0, 1000), nil)
		sout1, sout2 := string(alignment[0]), string(alignment[1])

		if sout1 != test.out1 || sout2 != test.out2 {
			t.Fatalf(
				`Alignment for: (sequence identitiy: %d)
%s
%s
%s
%s
resulted in
%s
%s
%s
%s
but should have been
%s
%s
%s
%s`,
				cablastp.SeqIdentity(alignment[0], alignment[1]),
				sep, test.seq1, test.seq2, sep,
				sep, sout1, sout2, sep,
				sep, test.out1, test.out2, sep)
		}
	}
}

func TestBiogoNeedlemanWunsch(t *testing.T) {
	type test struct {
		seq1, seq2 string
		out1, out2 string
	}

	tests := []test{
		{
			"ABCD",
			"ABCD",
			"ABCD",
			"ABCD",
		},
		{
			"GHIKLMNPQR",
			"GAAAHIKLMN",
			"---GHIKLMNPQR",
			"GAAAHIKLMN---",
		},
		{
			"GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
			"---GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
		},
		// { 
		// "ABCDEFGWXYZ", 
		// "ABCDEFMNPQRSTZABEGWXYZ", 
		// "ABCDEF-----------GWXYZ", 
		// "ABCDEFMNPQRSTZABEGWXYZ", 
		// }, 
	}
	aligner := &nw.Aligner{
		Matrix:  blosum.Matrix62,
		LookUp:  nwLookUpP,
		GapChar: '-',
	}
	sep := strings.Repeat("-", 45)
	for _, test := range tests {
		alignment, err := aligner.Align(
			&seq.Seq{Seq: []byte(test.seq1)}, &seq.Seq{Seq: []byte(test.seq2)})
		if err != nil {
			t.Fatal(err)
		}
		sout1, sout2 := string(alignment[0].Seq), string(alignment[1].Seq)

		if sout1 != test.out1 || sout2 != test.out2 {
			t.Fatalf(
				`Alignment for: (sequence identitiy: %d)
%s
%s
%s
%s
resulted in
%s
%s
%s
%s
but should have been
%s
%s
%s
%s`,
				cablastp.SeqIdentity(alignment[0].Seq, alignment[1].Seq),
				sep, test.seq1, test.seq2, sep,
				sep, sout1, sout2, sep,
				sep, test.out1, test.out2, sep)
		}
	}
}

// flags: flagMatchKmerSize, flagUngappedWindowSize, flagSeqIdThreshold,
// flagGappedWindowSize.
func TestExtendMatch(t *testing.T) {
	flagMatchKmerSize = 3
	flagUngappedWindowSize = 10
	flagSeqIdThreshold = 50
	flagGappedWindowSize = 25

	type test struct {
		rseq, oseq   string
		rmseq, omseq string
	}
	tests := []test{
		{
			"ABCDEFGHIKLMNPQR",
			"ABCDEFGHIKLMNPQR",
			"ABCDEFGHIKLMNPQR",
			"ABCDEFGHIKLMNPQR",
		},
		{
			"ABCDEFGHIKLMNPQRSTVW",
			"ABCDEFGAAAHIKLMNPQRSTVW",
			"ABCDEFGHIKLMNPQRSTVW",
			"ABCDEFGAAAHIKLMNPQRSTVW",
		},
		{
			"ABCDEFGHIKLMNPQRSTVW",
			"ABCDEFGAAAHIKLMNPQRSTBBBBBBBBBBBBBBBBBBBVW",
			"ABCDEF",
			"ABCDEF",
		},
	}
	sep := strings.Repeat("-", 45)
	for _, test := range tests {
		refMatch, orgMatch := extendMatch(
			[]byte(test.rseq), []byte(test.oseq),
			make([]byte, 0, 1000), make([]byte, 0, 1000), nil)
		srefMatch, sorgMatch := string(refMatch), string(orgMatch)

		if srefMatch != test.rmseq || sorgMatch != test.omseq {
			t.Fatalf(
				`Extending a match for:
%s
%s
%s
%s
resulted in
%s
%s
%s
%s
but should have been
%s
%s
%s
%s`,
				sep, test.rseq, test.oseq, sep,
				sep, srefMatch, sorgMatch, sep,
				sep, test.rmseq, test.omseq, sep)
		}
	}
}

// flags: flagMatchKmerSize, flagUngappedWindowSize, flagSeqIdThreshold
func TestUngappedExtension(t *testing.T) {
	flagMatchKmerSize = 3
	flagUngappedWindowSize = 10
	flagSeqIdThreshold = 50

	type test struct {
		rseq, oseq string
		answer     int
	}
	tests := []test{
		{"A", "A", 0},
		{"AB", "AB", 0},
		{"ABC", "ABC", 3},
		{"ABCD", "ABCD", 3},
		{"ABCYEFG", "ABCZEFG", 3},
		{"ABCYEFGH", "ABCZEFGH", 8},
		{"ABCDEFGHIJKLMNOP", "ABCDEFGHIJKLMNOP", 15},
		{"ABCDEF", "ABC", 3},
		{"ABC", "ABCDEF", 3},
	}

	for _, test := range tests {
		tval := alignUngapped([]byte(test.rseq), []byte(test.oseq))
		if tval != test.answer {
			t.Fatalf("Ungapped extension on '%s' and '%s' should yield a "+
				"length of %d, but 'alignUngapped' returned %d.",
				test.rseq, test.oseq, test.answer, tval)
		}
	}
}
