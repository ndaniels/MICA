package main

import (
	"strings"
	"testing"

	"github.com/BurntSushi/cablastp"
)

func TestSkipLowComplexity(t *testing.T) {
	type test struct {
		seq                    string
		windowSize, regionSize int
		skipped                int
		leftover               string
	}
	tests := []test{
		{"ABCDDDDDDDDDDDDDDDDDDXYZ", 10, 5, 21, "XYZ"},
		{"DDDDDDABCDEF", 10, 5, 6, "ABCDEF"},
		{"DDDDDDABCDEFFFFFFFFFFXYZXYZ", 10, 5, 6, "ABCDEFFFFFFFFFFXYZXYZ"},
		{"ABCDEFFFFFFFFFFFFFFFXYZXYZ", 10, 5, 20, "XYZXYZ"},
	}
	for _, test := range tests {
		skipped := skipLowComplexity(
			[]byte(test.seq), test.windowSize, test.regionSize)
		if skipped != test.skipped {
			t.Fatalf("Skipping low complexity regions in '%s' with a "+
				"window size of %d and a region size of %d should have "+
				"skipped %d residues, but it actually skipped %d residues.",
				test.seq, test.windowSize, test.regionSize,
				test.skipped, skipped)
		}

		leftover := test.seq[skipped:]
		if leftover != test.leftover {
			t.Fatalf("Skipping low complexity regions in '%s' with a "+
				"window size of %d and a region size of %d skipped %d "+
				"residues and returned '%s' as leftovers but should have "+
				"returned '%s'.",
				test.seq, test.windowSize, test.regionSize,
				test.skipped, leftover, test.leftover)
		}
	}
}

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
			"PPPGHIKLMNPQR",
			"GAAAHIKLMN",
			"PPPGHIKLMNPQR",
			"---GAAAHIKLMN",
		},
		{
			"GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
			"---GHIKLMNPQRSTVW",
			"GAAAHIKLMNPQRSTVW",
		},
		{
			"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
			"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
			"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
			"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
		},
		{
			"NNNNNNNN",
			"NNNNNNNN",
			"NNNNNNNN",
			"NNNNNNNN",
		},
		{
			"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
			"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
			"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
			"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
		},
		// {
		// "ABCDEFGWXYZ",
		// "ABCDEFMNPQRSTZABEGWXYZ",
		// "ABCDEF-----------GWXYZ",
		// "ABCDEFMNPQRSTZABEGWXYZ",
		// },
	}
	sep := strings.Repeat("-", 45)
	mem := newMemory()
	for _, test := range tests {
		alignment := nwAlign([]byte(test.seq1), []byte(test.seq2), mem)
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

func TestExtendMatch(t *testing.T) {
	flagMatchKmerSize := 3
	flagUngappedWindowSize := 10
	flagExtSeqIdThreshold := 50
	flagGappedWindowSize := 25

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
	mem := newMemory()
	for _, test := range tests {
		corMatch, orgMatch := extendMatch(
			[]byte(test.rseq), []byte(test.oseq),
			flagGappedWindowSize, flagUngappedWindowSize,
			flagMatchKmerSize, flagExtSeqIdThreshold,
			mem)
		scorMatch, sorgMatch := string(corMatch), string(orgMatch)

		if scorMatch != test.rmseq || sorgMatch != test.omseq {
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
				sep, scorMatch, sorgMatch, sep,
				sep, test.rmseq, test.omseq, sep)
		}
	}
}

func TestUngappedExtension(t *testing.T) {
	flagMatchKmerSize := 3
	flagUngappedWindowSize := 10
	flagExtSeqIdThreshold := 50

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
		{"ABCDEFGHIKLMNPQR", "ABCDEFGHIKLMNPQR", 15},
	}

	for _, test := range tests {
		tval := alignUngapped(
			[]byte(test.rseq), []byte(test.oseq),
			flagUngappedWindowSize, flagMatchKmerSize, flagExtSeqIdThreshold)
		if tval != test.answer {
			t.Fatalf("Ungapped extension on '%s' and '%s' should yield a "+
				"length of %d, but 'alignUngapped' returned %d.",
				test.rseq, test.oseq, test.answer, tval)
		}
	}
}
