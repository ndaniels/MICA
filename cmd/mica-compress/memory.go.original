package main

const (
	memSeqSize       = 10000
	dynamicTableSize = memSeqSize * memSeqSize
	numSeeds         = 100
)

// memory is a goroutine-specific memory arena, specifically used in each
// compression worker goroutine. Its purpose is to reduce the amount of
// memory allocation in hot-spots: sequence alignment and seed lookup.
type memory struct {
	table    []int
	ref, org []byte
	seeds    [][2]uint
}

func newMemory() *memory {
	return &memory{
		table: make([]int, memSeqSize*memSeqSize),
		ref:   make([]byte, memSeqSize),
		org:   make([]byte, memSeqSize),
		seeds: make([][2]uint, 0, numSeeds),
	}
}
