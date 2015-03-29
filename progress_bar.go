package neutronium

import (
	"sync/atomic"
)

type ProgressBar struct {
	Label   string
	Total   uint64
	Current uint64
}

func (bar *ProgressBar) Increment() {
	atomic.AddUint64(&bar.Current, 1)
}

func (bar *ProgressBar) ClearAndDisplay() {
	Vprint("\r")
	barWidth := uint64(80 - len(bar.Label))
	ticks := (barWidth * bar.Current) / bar.Total
	Vprintf("%s [", bar.Label)
	for i := uint64(0); i < ticks; i++ {
		Vprint("=")
	}
	for i := uint64(0); i < (barWidth - ticks); i++ {
		Vprint(" ")
	}
	Vprint("] ")
	Vprintf("%d / %d", bar.Current, bar.Total)
}
