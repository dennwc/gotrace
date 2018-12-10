package gotrace

func signf(v float64) float64 {
	if v > 0 {
		return +1
	} else if v < 0 {
		return -1
	} else {
		return 0
	}
}

func sign(v int) int {
	if v > 0 {
		return +1
	} else if v < 0 {
		return -1
	} else {
		return 0
	}
}

func abs(v int) int {
	if v >= 0 {
		return v
	}
	return -v
}

func fabs(v float64) float64 {
	if v >= 0 {
		return v
	}
	return -v
}

func mod(a, n int) int {
	if a >= n {
		return a % n
	} else if a >= 0 {
		return a
	}
	return n - 1 - (-1-a)%n
}

func floordiv(a, n int) int {
	if a >= 0 {
		return a / n
	}
	return -1 - (-1-a)/n
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// range over the straight line segment [a,b] when lambda ranges over [0,1]
func interval(lambda float64, a, b Point) (res Point) {
	res.X = a.X + lambda*(b.X-a.X)
	res.Y = a.Y + lambda*(b.Y-a.Y)
	return
}
