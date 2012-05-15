package main

import "math"

func pow(x, y int) int {
	return int(math.Pow(float64(x), float64(y)))
}

