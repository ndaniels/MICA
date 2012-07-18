package main

import "math"

func pow(x, y int) int {
	return int(math.Pow(float64(x), float64(y)))
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
