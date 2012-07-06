package main

type match struct {
	rseq *referenceSeq
	link *linkEntry
}

// Less implements an ordinal relationship between two link entries.
// Currently, one linkEntry is "better" than the other if it corresponds to a
// bigger part of a reference sequence in the compressed database.
func (m1 match) Less(m2 match) bool {
	le1, le2 := m1.link, m2.link
	return (le1.refEndRes - le1.refStartRes) < (le2.refEndRes - le2.refStartRes)
}
