package cablastp

import (
	"fmt"
	"log"
	"strconv"
	"strings"
)

const (
	ModSubstitution = iota
	ModDeletion
	ModInsertion
)

type EditScript struct {
	Mods []*Mod
}

func NewEditScript(alignment [2][]byte) *EditScript {
	return newEditScript(alignment[0], alignment[1])
}

func NewEditScriptParse(editScript string) (*EditScript, error) {
	mods := make([]*Mod, 0, 15)
	var mod *Mod = nil

	for i := 0; i < len(editScript); i++ {
		b := editScript[i]

		switch b {
		case 's':
			fallthrough
		case 'i':
			fallthrough
		case 'd':
			if mod != nil {
				mod.End = mod.Start + len(mod.Residues)
				mods = append(mods, mod)
			}
			newMod := NewMod(byteToModKind(b))

			// Consume numbers until we hit a non-number.
			digits := make([]byte, 0, 3)
			for j := i + 1; j < len(editScript); j++ {
				if editScript[j] >= '0' && editScript[j] <= '9' {
					digits = append(digits, editScript[j])
				} else {
					i = j - 1
					break
				}
			}

			// If we didn't find any numbers, we have a syntax error.
			if len(digits) == 0 {
				return nil, fmt.Errorf("Expected an offset number after '%c' "+
					"in column %d of '%s'.", b, i, editScript)
			}

			// If we can't parse the number as an integer, syntax error!
			num64, err := strconv.ParseInt(string(digits), 10, 32)
			if err != nil {
				return nil, fmt.Errorf("Expected an offset number after '%c' "+
					"in column %d of '%s', but got '%s' instead.",
					b, i, editScript, string(digits))
			}

			newMod.Start = int(num64)
			if mod != nil {
				newMod.Start += mod.Start
			}
			mod = newMod
		default:
			if b >= '0' && b <= '9' {
				return nil, fmt.Errorf("Expected a residue at column %d in "+
					"'%s', but got a number '%c' instead.", i, editScript, b)
			}
			if mod == nil {
				return nil, fmt.Errorf("Expected 's', 'i' or 'd' but got '%c' "+
					"at column %d in '%s'.", b, i, editScript)
			}

			mod.Residues = append(mod.Residues, b)
		}
	}

	// One last modification?
	if mod != nil {
		mod.End = mod.Start + len(mod.Residues)
		mods = append(mods, mod)
	}

	return &EditScript{
		Mods: mods,
	}, nil
}

func newEditScript(fromSeq, toSeq []byte) *EditScript {
	if len(fromSeq) != len(toSeq) {
		log.Panicf("A new edit script can only be generated with two "+
			"sequences of equal length. Lengths of %d and %d were provided.",
			len(fromSeq), len(toSeq))
	}

	// The set of all modifications between fromSeq and toSeq, expressed
	// in terms of substitutions, insertions and deletions.
	mods := make([]*Mod, 0, 15)

	// mod corresponds to the current modification that is taking
	// place. It is either nil (no modification), a substitution, a deletion
	// or an insertion.
	var mod *Mod = nil

	// fromIndex corresponds to the true index of the 'fromSeq'.
	// This is necessary because an alignment puts '-' characters in 'fromSeq',
	// which we don't want to count towards fromSeq index information.
	fromIndex := 0

	for i := 0; i < len(fromSeq); i++ {
		from, to := fromSeq[i], toSeq[i]

		newModKind := -1
		switch {
		case from == to:
			newModKind = -1
		case from == '-':
			newModKind = ModInsertion
		case to == '-':
			newModKind = ModDeletion
		case from != to:
			newModKind = ModSubstitution
		default:
			log.Panicf("BUG: from: %c, to: %c", from, to)
		}

		if mod != nil {
			if mod.Kind == newModKind {
				mod.AddResidue(to)
			} else {
				mod.End = fromIndex
				mods = append(mods, mod)

				if newModKind == -1 {
					mod = nil
				} else {
					mod = NewMod(newModKind)
					mod.Start = fromIndex
					mod.End = fromIndex // never changes if mod is ModInsertion
					mod.AddResidue(to)
				}
			}
		} else if newModKind != -1 {
			mod = NewMod(newModKind)
			mod.Start = fromIndex
			mod.End = fromIndex // never changes if mod is ModInsertion
			mod.AddResidue(to)
		}

		if from != '-' {
			fromIndex++
		}
	}

	return &EditScript{
		Mods: mods,
	}
}

func (diff *EditScript) Apply(fromSeq []byte) []byte {
	toSeq := make([]byte, 0, len(fromSeq))
	lastEnd := 0
	for _, mod := range diff.Mods {
		toSeq = append(toSeq, fromSeq[lastEnd:mod.Start]...)
		toSeq = append(toSeq, mod.Residues...)
		lastEnd = mod.End
	}
	toSeq = append(toSeq, fromSeq[lastEnd:len(fromSeq)]...)
	return toSeq
}

func (diff *EditScript) String() string {
	mods := make([]string, len(diff.Mods))
	lastDist := 0
	for i, m := range diff.Mods {
		dist := m.Start - lastDist
		lastDist = m.Start

		switch m.Kind {
		case ModSubstitution:
			mods[i] = fmt.Sprintf("s%d%s", dist,
				strings.ToUpper(string(m.Residues)))
		case ModInsertion:
			mods[i] = fmt.Sprintf("i%d%s", dist,
				strings.ToUpper(string(m.Residues)))
		case ModDeletion:
			mods[i] = fmt.Sprintf("d%d%s", dist,
				strings.Repeat("-", m.End-m.Start))
		default:
			log.Panicf("Invalid kind '%d' for an EditScript modification.",
				m.Kind)
		}
	}
	return strings.Join(mods, "")
}

type Mod struct {
	// Either a substitution, deletion or an insertion
	Kind int

	// Indices into the 'from' sequence where the modification starts and ends.
	// These are equivalent for insertions.
	Start, End int

	// The new residues.
	// When Mod is a deletion, this is empty. In which case, the number of 
	// residues to delete equals (End - Start).
	Residues []byte
}

func NewMod(kind int) *Mod {
	return &Mod{
		Kind:     kind,
		Start:    0,
		End:      0,
		Residues: make([]byte, 0, 20),
	}
}

func (m *Mod) AddResidue(residue byte) {
	if m.Kind != ModDeletion {
		m.Residues = append(m.Residues, residue)
	}
}

func (m *Mod) String() string {
	switch m.Kind {
	case ModSubstitution:
		return fmt.Sprintf("s(%d,%d)%s", m.Start, m.End, string(m.Residues))
	case ModInsertion:
		return fmt.Sprintf("i(%d,%d)%s", m.Start, m.End, string(m.Residues))
	case ModDeletion:
		return fmt.Sprintf("d(%d,%d)%s", m.Start, m.End, string(m.Residues))
	}
	log.Panicf("Invalid kind '%d' for an EditScript modification.", m.Kind)
	panic("unreachable")
}

func byteToModKind(char byte) int {
	switch char {
	case 's':
		return ModSubstitution
	case 'i':
		return ModInsertion
	case 'd':
		return ModDeletion
	}
	panic("unreachable")
}
