package neutronium

func reduce1(char byte) byte {
	switch char {
	case 'F', 'W', 'Y':
		return 'A'
	case 'C', 'I', 'L', 'M', 'V', 'J':
		return 'C'
	case 'A', 'G', 'P', 'S', 'T':
		return 'G'
	case 'D', 'E', 'N', 'Q', 'K', 'R', 'H', 'B', 'Z':
		return 'T'
	case 'X':
		return 'N'
	}
	panic("bad letter, could not reduce")
}

func Reduce(seq []byte) []byte {
	dest := make([]byte, len(seq))
	for i, residue := range seq {
		dest[i] = reduce1(residue)
	}
	return dest
}
