package cablastp

import (
	"bytes"
	"strings"
	"testing"
)

func TestDBConfIO(t *testing.T) {
	dbConf := DefaultDBConf
	buf := new(bytes.Buffer)

	dbConf.Write(buf)
	dbConfTest, err := LoadDBConf(buf)
	if err != nil {
		t.Fatal(err)
	}
	if dbConf != dbConfTest {
		t.Fatalf("%v != %v", dbConf, dbConfTest)
	}
}

func TestEditScripts(t *testing.T) {
	type test struct {
		fromSeq, toSeq               string
		alignedFromSeq, alignedToSeq string
		script                       string
	}

	tests := []test{
		{
			// Yes, this is a nucleotide sequence. But that doesn't matter.
			// This is taken straight from the CaBLAST paper.
			"GTTCACTTATGTATTCATATGATTTTGGCAA",
			"GTTCACGTGTATATTTATATAATTTTGGCAA",
			"GTTCACTTATGTATTC--ATATGATTTTGGCAA", // in alignment
			"GTTCACG--TGTATATTTATATAATTTTGGCAA", // in alignment
			"s6Gd1--s7ATi2TTs4A",
		},
	}
	sep := strings.Repeat("-", 45)
	for _, test := range tests {
		// Check that the edit script generated by the alignment of 'from'
		// and 'to' is correct.
		tval := newEditScript(
			[]byte(test.alignedFromSeq), []byte(test.alignedToSeq))

		t.Log("")
		for _, mod := range tval.Mods {
			t.Log(mod.String())
		}

		if tval.String() != test.script {
			t.Fatalf(
				`The EditScript generated by
%s
%s
%s
%s
resulted in
%s
%s
%s
but should be
%s
%s
%s`,
				sep, test.alignedFromSeq, test.alignedToSeq, sep,
				sep, tval.String(), sep,
				sep, test.script, sep)
		}

		// Check that parsing the script and converting it back to a string
		// yields the same script.
		checkScript, err := NewEditScriptParse(test.script)
		if err != nil {
			t.Fatalf("Parsing '%s' caused: %s", test.script, err)
		}
		if checkScript.String() != test.script {
			t.Fatalf(
				`Parsing
%s
%s
%s
and converting it back to a string resulted in
%s
%s
%s
but should be
%s
%s
%s`,
				sep, test.script, sep,
				sep, checkScript.String(), sep,
				sep, test.script, sep)
		}

		// Finally, check that re-applying test.script to fromSeq
		// yields toSeq exactly.
		toSeq := string(tval.Apply([]byte(test.fromSeq)))
		if toSeq != test.toSeq {
			t.Fatalf("Applying the edit script '%s' to the sequence '%s' "+
				"resulted in\n%s\n%s\n%s\nbut should be\n%s\n%s\n%s",
				test.script, test.fromSeq,
				sep, toSeq, sep,
				sep, test.toSeq, sep)
		}
	}
}
