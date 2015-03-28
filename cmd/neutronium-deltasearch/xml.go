package main

import (
	"encoding/xml"
)

type blast struct {
	XMLName xml.Name `xml:"BlastOutput"`
	Hits    []hit    `xml:"BlastOutput_iterations>Iteration>Iteration_hits>Hit"`
}

type hit struct {
	XMLName   xml.Name `xml:"Hit"`
	Num       int      `xml:"Hit_num"`
	Accession int      `xml:"Hit_accession"`
	Hsps      []hsp    `xml:"Hit_hsps>Hsp"`
}

type hsp struct {
	XMLName   xml.Name `xml:"Hsp"`
	Num       int      `xml:"Hsp_num"`
	Evalue    float64  `xml:"Hsp_evalue"`
	QueryFrom int      `xml:"Hsp_query-from"`
	QueryTo   int      `xml:"Hsp_query-to"`
	HitFrom   int      `xml:"Hsp_hit-from"`
	HitTo     int      `xml:"Hsp_hit-to"`
}
