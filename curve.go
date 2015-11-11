package gotrace

type segment struct {
	Segment

	vertex Point
	alpha  float64
	alpha0 float64
	beta   float64
}

func newPrivCurve(n int) privCurve {
	return privCurve{segm: make([]segment, n)}
}

type privCurve struct {
	alphacurve bool
	segm       []segment
}

func (c *privCurve) ToCurve() []Segment {
	out := make([]Segment, len(c.segm))
	for i, s := range c.segm {
		out[i] = s.Segment
	}
	return out
}
