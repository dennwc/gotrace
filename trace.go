package gotrace

import (
	"fmt"
	"math"
)

const (
	cINFTY  = 10000000        // it suffices that this is longer than any path; it need not be really infinite
	cCOS179 = -0.999847695156 // the cosine of 179 degrees
)

// return a direction that is 90 degrees counterclockwise from p2-p0,
// but then restricted to one of the major wind directions (n, nw, w, etc)
func dorthInfty(p0, p2 Point) Point {
	return Point{signf(p2.X - p0.X), -signf(p2.Y - p0.Y)}
}

// return (p1-p0)x(p2-p0), the area of the parallelogram
func dpara(p0, p1, p2 Point) float64 {
	x1 := p1.X - p0.X
	y1 := p1.Y - p0.Y
	x2 := p2.X - p0.X
	y2 := p2.Y - p0.Y
	return x1*y2 - x2*y1
}

// ddenom/dpara have the property that the square of radius 1 centered
// at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2)
func ddenom(p0, p2 Point) float64 {
	r := dorthInfty(p0, p2)
	return r.Y*(p2.X-p0.X) - r.X*(p2.Y-p0.Y)
}

// return 1 if a <= b < c < a, in a cyclic sense (mod n)
func cyclic(a, b, c int) bool {
	if a <= c {
		return a <= b && b < c
	}
	return a <= b || b < c
}

type sums struct {
	x, y       int
	x2, xy, y2 int
}

type point struct {
	X, Y int
}

// the path structure is filled in with information about a given path
// as it is accumulated and passed through the different stages of the
// Potrace algorithm. Backends only need to read the fcurve and fm
// fields of this data structure, but debugging backends may read
// other fields.
type privPath struct {
	Pt  []point // path as extracted from bitmap
	Lon []int   // (i,lon[i]) = longest straight line from i

	Orig point  // origin for sums
	Sums []sums // sums[len+1]: cache for fast summing

	Po []int // optimal polygon

	Curve  privCurve  // curve[m]: array of curve elements
	OCurve privCurve  // ocurve[om]: array of curve elements
	FCurve *privCurve // final curve: this points to either curve or ocurve
}

// determine the center and slope of the line i..j. Assume i<j. Needs
// "sum" components of p to be set.
func pointslope(pp *privPath, i, j int) (ctr, dir Point) {
	// assume i<j

	n := len(pp.Pt)
	sums := pp.Sums
	r := 0 // rotations from i to j

	for j >= n {
		j -= n
		r++
	}
	for i >= n {
		i -= n
		r--
	}
	for j < 0 {
		j += n
		r--
	}
	for i < 0 {
		i += n
		r++
	}

	x := float64(sums[j+1].x - sums[i].x + r*sums[n].x)
	y := float64(sums[j+1].y - sums[i].y + r*sums[n].y)
	x2 := float64(sums[j+1].x2 - sums[i].x2 + r*sums[n].x2)
	xy := float64(sums[j+1].xy - sums[i].xy + r*sums[n].xy)
	y2 := float64(sums[j+1].y2 - sums[i].y2 + r*sums[n].y2)
	k := float64(j + 1 - i + r*n)

	ctr.X = x / k
	ctr.Y = y / k

	a := (x2 - x*x/k) / k
	b := (xy - x*y/k) / k
	c := (y2 - y*y/k) / k

	lambda2 := (a + c + math.Sqrt((a-c)*(a-c)+4*b*b)) / 2 // larger e.value

	// now find e.vector for lambda2
	a -= lambda2
	c -= lambda2

	var l float64
	if fabs(a) >= fabs(c) {
		l = math.Sqrt(a*a + b*b)
		if l != 0 {
			dir.X = -b / l
			dir.Y = a / l
		}
	} else {
		l = math.Sqrt(c*c + b*b)
		if l != 0 {
			dir.X = -c / l
			dir.Y = b / l
		}
	}
	if l == 0 {
		dir.X, dir.Y = 0, 0 // sometimes this can happen when k=4: the two eigenvalues coincide
	}
	return
}

// calculate p1 x p2
func xprod(p1, p2 point) int {
	return int(p1.X*p2.Y - p1.Y*p2.X)
}

// calculate (p1-p0)x(p3-p2)
func cprod(p0, p1, p2, p3 Point) float64 {
	x1 := p1.X - p0.X
	y1 := p1.Y - p0.Y
	x2 := p3.X - p2.X
	y2 := p3.Y - p2.Y
	return x1*y2 - x2*y1
}

// calculate (p1-p0)*(p2-p0)
func iprod(p0, p1, p2 Point) float64 {
	x1 := p1.X - p0.X
	y1 := p1.Y - p0.Y
	x2 := p2.X - p0.X
	y2 := p2.Y - p0.Y
	return x1*x2 + y1*y2
}

// calculate (p1-p0)*(p3-p2)
func iprod1(p0, p1, p2, p3 Point) float64 {
	x1 := p1.X - p0.X
	y1 := p1.Y - p0.Y
	x2 := p3.X - p2.X
	y2 := p3.Y - p2.Y
	return x1*x2 + y1*y2
}

// calculate distance between two points
func ddist(p, q Point) float64 {
	x, y := p.X-q.X, p.Y-q.Y
	return math.Sqrt(x*x + y*y)
}

// calculate point of a bezier curve
func bezier(t float64, p0, p1, p2, p3 Point) (res Point) {
	s := 1 - t
	res.X = s*s*s*p0.X + 3*(s*s*t)*p1.X + 3*(t*t*s)*p2.X + t*t*t*p3.X
	res.Y = s*s*s*p0.Y + 3*(s*s*t)*p1.Y + 3*(t*t*s)*p2.Y + t*t*t*p3.Y
	return
}

// calculate the point t in [0..1] on the (convex) bezier curve (p0,p1,p2,p3)
// which is tangent to q1-q0. Return -1.0 if there is no solution in [0..1].
func tangent(p0, p1, p2, p3, q0, q1 Point) float64 {
	A := cprod(p0, p1, q0, q1)
	B := cprod(p1, p2, q0, q1)
	C := cprod(p2, p3, q0, q1)

	a := A - 2*B + C
	b := -2*A + 2*B
	c := A

	d := b*b - 4*a*c

	if a == 0 || d < 0 {
		return -1 // TODO: replace with error?
	}
	s := math.Sqrt(d)

	r1 := (-b + s) / (2 * a)
	r2 := (-b - s) / (2 * a)

	if r1 >= 0 && r1 <= 1 {
		return r1
	} else if r2 >= 0 && r2 <= 1 {
		return r2
	} else {
		return -1
	}
}

// the type of (affine) quadratic forms, represented as symmetric 3x3
// matrices.  The value of the quadratic form at a vector (x,y) is v^t
// Q v, where v = (x,y,1)^t.
type quadForm [3][3]float64

// Apply quadratic form Q to vector w = (w.x,w.y)
func quadform(Q quadForm, w Point) (sum float64) {
	v := [3]float64{w.X, w.Y, 1}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			sum += v[i] * Q[i][j] * v[j]
		}
	}
	return sum
}

// calcSums is a preparation: fill in the sum* fields of a path (used for later rapid summing).
func (pp *privPath) calcSums() {
	pp.Sums = make([]sums, len(pp.Pt)+1)
	pp.Orig = pp.Pt[0]
	for i, p := range pp.Pt {
		x := p.X - pp.Orig.X
		y := p.Y - pp.Orig.Y
		pp.Sums[i+1].x = pp.Sums[i].x + x
		pp.Sums[i+1].y = pp.Sums[i].y + y
		pp.Sums[i+1].x2 = pp.Sums[i].x2 + x*x
		pp.Sums[i+1].xy = pp.Sums[i].xy + x*y
		pp.Sums[i+1].y2 = pp.Sums[i].y2 + y*y
	}
}

/* Stage 1: determine the straight subpaths (Sec. 2.2.1). Fill in the
   "lon" component of a path object (based on pt/len).	For each i,
   lon[i] is the furthest index such that a straight line can be drawn
   from i to lon[i]. Return 1 on error with errno set, else 0. */

/* this algorithm depends on the fact that the existence of straight
   subpaths is a triplewise property. I.e., there exists a straight
   line through squares i0,...,in iff there exists a straight line
   through i,j,k, for all i0<=i<j<k<=in. (Proof?) */

/* this implementation of calcLon is O(n^2). It replaces an older
   O(n^3) version. A "constraint" means that future points must
   satisfy xprod(constraint[0], cur) >= 0 and xprod(constraint[1],
   cur) <= 0. */

/* Remark for Potrace 1.1: the current implementation of calcLon is
   more complex than the implementation found in Potrace 1.0, but it
   is considerably faster. The introduction of the "nc" data structure
   means that we only have to test the constraints for "corner"
   points. On a typical input file, this speeds up the calcLon
   function by a factor of 31.2, thereby decreasing its time share
   within the overall Potrace algorithm from 72.6% to 7.82%, and
   speeding up the overall algorithm by a factor of 3.36. On another
   input file, calcLon was sped up by a factor of 6.7, decreasing its
   time share from 51.4% to 13.61%, and speeding up the overall
   algorithm by a factor of 1.78. In any case, the savings are
   substantial. */

func (pp *privPath) calcLon() {
	pt := pp.Pt
	n := len(pt)
	var (
		pivk = make([]int, n)
		nc   = make([]int, n) // next corner
	)

	// 	 initialize the nc data structure. Point from each point to the
	// 	 furthest future point to which it is connected by a vertical or
	// 	 horizontal segment. We take advantage of the fact that there is
	// 	 always a direction change at 0 (due to the path decomposition
	// 	 algorithm). But even if this were not so, there is no harm, as
	// 	 in practice, correctness does not depend on the word "furthest"
	// 	 above.
	k := 0
	for i := n - 1; i >= 0; i-- {
		if pt[i].X != pt[k].X && pt[i].Y != pt[k].Y {
			k = i + 1 // necessarily i<n-1 in this case
		}
		nc[i] = k
	}

	pp.Lon = make([]int, n)

	// determine pivot points: for each i, let pivk[i] be the furthest k
	//  such that all j with i<j<k lie on a line connecting i,k.

	var (
		ct         [4]int
		constraint [2]point
		cur, off   point
		dk         point
		j          int

		a, b, c, d int
	)
	for i := n - 1; i >= 0; i-- {
		ct[0], ct[1], ct[2], ct[3] = 0, 0, 0, 0

		// keep track of "directions" that have occurred
		dir := int((3 + 3*(pt[mod(i+1, n)].X-pt[i].X) + (pt[mod(i+1, n)].Y - pt[i].Y)) / 2)
		ct[dir]++

		constraint[0].X = 0
		constraint[0].Y = 0
		constraint[1].X = 0
		constraint[1].Y = 0

		// find the next k such that no straight line from i to k
		k = nc[i]
		k1 := i
		for {
			dir = (3 + 3*int(sign(pt[k].X-pt[k1].X)) + int(sign(pt[k].Y-pt[k1].Y))) / 2
			ct[dir]++

			// if all four "directions" have occurred, cut this path
			if ct[0] != 0 && ct[1] != 0 && ct[2] != 0 && ct[3] != 0 {
				pivk[i] = k1
				goto foundk
			}

			cur.X = pt[k].X - pt[i].X
			cur.Y = pt[k].Y - pt[i].Y

			// see if current constraint is violated
			if xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0 {
				goto constraint_viol
			}

			// else, update constraint
			if abs(cur.X) <= 1 && abs(cur.Y) <= 1 {
				// no constraint
			} else {
				if cur.Y >= 0 && (cur.Y > 0 || cur.X < 0) {
					off.X = cur.X + 1
				} else {
					off.X = cur.X - 1
				}
				if cur.X <= 0 && (cur.X < 0 || cur.Y < 0) {
					off.Y = cur.Y + 1
				} else {
					off.Y = cur.Y - 1
				}
				if xprod(constraint[0], off) >= 0 {
					constraint[0] = off
				}
				if cur.Y <= 0 && (cur.Y < 0 || cur.X < 0) {
					off.X = cur.X + 1
				} else {
					off.X = cur.X - 1
				}
				if cur.X >= 0 && (cur.X > 0 || cur.Y < 0) {
					off.Y = cur.Y + 1
				} else {
					off.Y = cur.Y - 1
				}
				if xprod(constraint[1], off) <= 0 {
					constraint[1] = off
				}
			}
			k1 = k
			k = nc[k1]
			if !cyclic(k, i, k1) {
				break
			}
		}
	constraint_viol:
		//    k1 was the last "corner" satisfying the current constraint, and
		//    k is the first one violating it. We now need to find the last
		//    point along k1..k which satisfied the constraint.
		dk.X = sign(pt[k].X - pt[k1].X)
		dk.Y = sign(pt[k].Y - pt[k1].Y)
		cur.X = pt[k1].X - pt[i].X
		cur.Y = pt[k1].Y - pt[i].Y
		//    find largest integer j such that xprod(constraint[0], cur+j*dk)
		//    >= 0 and xprod(constraint[1], cur+j*dk) <= 0. Use bilinearity
		//    of xprod.
		a = xprod(constraint[0], cur)
		b = xprod(constraint[0], dk)
		c = xprod(constraint[1], cur)
		d = xprod(constraint[1], dk)
		/* find largest integer j such that a+j*b>=0 and c+j*d<=0. This
		   can be solved with integer arithmetic. */
		j = cINFTY
		if b < 0 {
			j = floordiv(a, -b)
		}
		if d > 0 {
			j = min(j, floordiv(-c, d))
		}
		pivk[i] = mod(k1+j, n)
	foundk:
	} // for i

	// clean up: for each i, let lon[i] be the largest k such that for
	//  all i' with i<=i'<k, i'<k<=pivk[i'].

	j = pivk[n-1]
	pp.Lon[n-1] = j
	for i := n - 2; i >= 0; i-- {
		if cyclic(i+1, pivk[i], j) {
			j = pivk[i]
		}
		pp.Lon[i] = j
	}

	for i := n - 1; cyclic(mod(i+1, n), j, pp.Lon[i]); i-- {
		pp.Lon[i] = j
	}
}

// Stage 2: calculate the optimal polygon (Sec. 2.2.2-2.2.4).

// Auxiliary function: calculate the penalty of an edge from i to j in
// the given path. This needs the "lon" and "sum*" data.
func penalty3(pp *privPath, i, j int) float64 {
	r := 0 // rotations from i to j
	sums, pt := pp.Sums, pp.Pt
	n := len(pt)

	if j >= n {
		j -= n
		r = 1
	}

	var (
		x, y, x2, xy, y2, k float64
	)

	// critical inner loop: the "if" gives a 4.6 percent speedup
	if r == 0 {
		x = float64(sums[j+1].x - sums[i].x)
		y = float64(sums[j+1].y - sums[i].y)
		x2 = float64(sums[j+1].x2 - sums[i].x2)
		xy = float64(sums[j+1].xy - sums[i].xy)
		y2 = float64(sums[j+1].y2 - sums[i].y2)
		k = float64(j + 1 - i)
	} else {
		x = float64(sums[j+1].x - sums[i].x + sums[n].x)
		y = float64(sums[j+1].y - sums[i].y + sums[n].y)
		x2 = float64(sums[j+1].x2 - sums[i].x2 + sums[n].x2)
		xy = float64(sums[j+1].xy - sums[i].xy + sums[n].xy)
		y2 = float64(sums[j+1].y2 - sums[i].y2 + sums[n].y2)
		k = float64(j + 1 - i + n)
	}

	px := float64(pt[i].X+pt[j].X)/2.0 - float64(pt[0].X)
	py := float64(pt[i].Y+pt[j].Y)/2.0 - float64(pt[0].Y)
	ey := float64(pt[j].X - pt[i].X)
	ex := -float64(pt[j].Y - pt[i].Y)

	a := ((x2-2*x*px)/k + px*px)
	b := ((xy-x*py-y*px)/k + px*py)
	c := ((y2-2*y*py)/k + py*py)

	s := ex*ex*a + 2*ex*ey*b + ey*ey*c

	return math.Sqrt(s)
}

// find the optimal polygon. Fill in the po component.
// Non-cyclic version: assumes i=0 is in the polygon. Fixme: implement cyclic version.
func (pp *privPath) bestpolygon() {
	//int i, j, m, k;
	n := len(pp.Pt)
	var (
		pen     = make([]float64, n+1) // penalty vector
		prev    = make([]int, n+1)     // best path pointer vector
		clip0   = make([]int, n)       // longest segment pointer, non-cyclic
		clip1   = make([]int, n+1)     // backwards segment pointer, non-cyclic
		seg0    = make([]int, n+1)     // forward segment bounds, m<=n
		seg1    = make([]int, n+1)     // backward segment bounds, m<=n
		thispen float64
		best    float64
		c       int
	)

	// calculate clipped paths
	for i := 0; i < n; i++ {
		c = mod(pp.Lon[mod(i-1, n)]-1, n)
		if c == i {
			c = mod(i+1, n)
		}
		if c < i {
			clip0[i] = n
		} else {
			clip0[i] = c
		}
	}

	// calculate backwards path clipping, non-cyclic. j <= clip0[i] iff clip1[j] <= i, for i,j=0..n.
	j := 1
	for i := 0; i < n; i++ {
		for j <= clip0[i] {
			clip1[j] = i
			j++
		}
	}

	// calculate seg0[j] = longest path from 0 with j segments
	i := 0
	for j = 0; i < n; j++ {
		seg0[j] = i
		i = clip0[i]
	}
	seg0[j] = n
	m := j

	// calculate seg1[j] = longest path to n with m-j segments
	i = n
	for j = m; j > 0; j-- {
		seg1[j] = i
		i = clip1[i]
	}
	seg1[0] = 0

	// now find the shortest path with m segments, based on penalty3
	// note: the outer 2 loops jointly have at most n iterations, thus
	//  the worst-case behavior here is quadratic. In practice, it is
	//  close to linear since the inner loop tends to be short.
	pen[0] = 0
	for j = 1; j <= m; j++ {
		for i = seg1[j]; i <= seg0[j]; i++ {
			best = -1
			for k := seg0[j-1]; k >= clip1[i]; k-- {
				thispen = penalty3(pp, k, i) + pen[k]
				if best < 0 || thispen < best {
					prev[i] = k
					best = thispen
				}
			}
			pen[i] = best
		}
	}

	pp.Po = make([]int, m)

	// read off shortest path
	for i, j = n, m-1; i > 0; j-- {
		i = prev[i]
		pp.Po[j] = i
	}
}

// Stage 3: vertex adjustment (Sec. 2.3.1).

// adjustVertices adjusts vertices of optimal polygon: calculate the intersection of
// the two "optimal" line segments, then move it into the unit square
// if it lies outside. Return 1 with errno set on error; 0 on success.
func (pp *privPath) adjustVertices() {
	po := pp.Po
	m := len(po)
	pt := pp.Pt
	n := len(pt)
	x0, y0 := pp.Orig.X, pp.Orig.Y

	var (
		ctr = make([]Point, m)
		dir = make([]Point, m)
		q   = make([]quadForm, m)
		v   [3]float64
		s   Point
	)

	pp.Curve = newPrivCurve(m)

	// calculate "optimal" point-slope representation for each line segment
	for i := 0; i < m; i++ {
		j := po[mod(i+1, m)]
		j = mod(j-po[i], n) + po[i]
		ctr[i], dir[i] = pointslope(pp, po[i], j)
	}

	// represent each line segment as a singular quadratic form; the
	// distance of a point (x,y) from the line segment will be
	// (x,y,1)Q(x,y,1)^t, where Q=q[i].
	for i := 0; i < m; i++ {
		d := dir[i].X*dir[i].X + dir[i].Y*dir[i].Y
		if d == 0.0 {
			for j := 0; j < 3; j++ {
				for k := 0; k < 3; k++ {
					q[i][j][k] = 0
				}
			}
		} else {
			v[0] = dir[i].Y
			v[1] = -dir[i].X
			v[2] = -v[1]*ctr[i].Y - v[0]*ctr[i].X
			for l := 0; l < 3; l++ {
				for k := 0; k < 3; k++ {
					q[i][l][k] = v[l] * v[k] / d
				}
			}
		}
	}

	// now calculate the "intersections" of consecutive segments.
	// Instead of using the actual intersection, we find the point
	// within a given unit square which minimizes the square distance to
	// the two lines.
	for i := 0; i < m; i++ {
		var (
			Q           quadForm
			w           Point
			dx, dy, det float64
			min, cand   float64 // minimum and candidate for minimum of quad. form
			xmin, ymin  float64 // coordinates of minimum
			z           int
		)

		// let s be the vertex, in coordinates relative to x0/y0
		s.X = float64(pt[po[i]].X - x0)
		s.Y = float64(pt[po[i]].Y - y0)

		// intersect segments i-1 and i

		j := mod(i-1, m)

		// add quadratic forms
		for l := 0; l < 3; l++ {
			for k := 0; k < 3; k++ {
				Q[l][k] = q[j][l][k] + q[i][l][k]
			}
		}

		for {
			// minimize the quadratic form Q on the unit square
			// find intersection

			det = Q[0][0]*Q[1][1] - Q[0][1]*Q[1][0]
			if det != 0.0 {
				w.X = (-Q[0][2]*Q[1][1] + Q[1][2]*Q[0][1]) / det
				w.Y = (Q[0][2]*Q[1][0] - Q[1][2]*Q[0][0]) / det
				break
			}

			// matrix is singular - lines are parallel. Add another,
			// orthogonal axis, through the center of the unit square
			if Q[0][0] > Q[1][1] {
				v[0] = -Q[0][1]
				v[1] = Q[0][0]
			} else if Q[1][1] != 0 {
				v[0] = -Q[1][1]
				v[1] = Q[1][0]
			} else {
				v[0] = 1
				v[1] = 0
			}
			d := v[0]*v[0] + v[1]*v[1]
			v[2] = -v[1]*float64(s.Y) - v[0]*float64(s.X)
			for l := 0; l < 3; l++ {
				for k := 0; k < 3; k++ {
					Q[l][k] += v[l] * v[k] / d
				}
			}
		}
		dx = fabs(w.X - s.X)
		dy = fabs(w.Y - s.Y)
		if dx <= .5 && dy <= .5 {
			pp.Curve.segm[i].vertex.X = w.X + float64(x0)
			pp.Curve.segm[i].vertex.Y = w.Y + float64(y0)
			continue
		}

		// the minimum was not in the unit square; now minimize quadratic on boundary of square
		min = quadform(Q, s)
		xmin = s.X
		ymin = s.Y

		if Q[0][0] == 0.0 {
			goto fixx
		}
		for z = 0; z < 2; z++ { // value of the y-coordinate
			w.Y = float64(s.Y) - 0.5 + float64(z)
			w.X = -(Q[0][1]*w.Y + Q[0][2]) / Q[0][0]
			dx = fabs(w.X - s.X)
			cand = quadform(Q, w)
			if dx <= .5 && cand < min {
				min = cand
				xmin = w.X
				ymin = w.Y
			}
		}
	fixx:
		if Q[1][1] == 0.0 {
			goto corners
		}
		for z = 0; z < 2; z++ { // value of the x-coordinate
			w.X = s.X - 0.5 + float64(z)
			w.Y = -(Q[1][0]*w.X + Q[1][2]) / Q[1][1]
			dy = fabs(w.Y - s.Y)
			cand = quadform(Q, w)
			if dy <= .5 && cand < min {
				min = cand
				xmin = w.X
				ymin = w.Y
			}
		}
	corners:
		// check four corners
		for l := 0; l < 2; l++ {
			for k := 0; k < 2; k++ {
				w.X = s.X - 0.5 + float64(l)
				w.Y = s.Y - 0.5 + float64(k)
				cand = quadform(Q, w)
				if cand < min {
					min = cand
					xmin = w.X
					ymin = w.Y
				}
			}
		}

		pp.Curve.segm[i].vertex.X = xmin + float64(x0)
		pp.Curve.segm[i].vertex.Y = ymin + float64(y0)
		continue
	}
}

// Stage 4: smoothing and corner analysis (Sec. 2.3.3)

// reverse orientation of a path
func reverse(curve *privCurve) {
	m := len(curve.segm)
	for i, j := 0, m-1; i < j; i, j = i+1, j-1 {
		curve.segm[i].vertex, curve.segm[j].vertex = curve.segm[j].vertex, curve.segm[i].vertex
	}
}

func smooth(curve *privCurve, alphamax float64) {
	m := len(curve.segm)
	// examine each vertex and find its best fit
	for i := 0; i < m; i++ {
		j := mod(i+1, m)
		k := mod(i+2, m)
		p4 := interval(1/2.0, curve.segm[k].vertex, curve.segm[j].vertex)

		var alpha float64
		denom := ddenom(curve.segm[i].vertex, curve.segm[k].vertex)
		if denom != 0.0 {
			dd := dpara(curve.segm[i].vertex, curve.segm[j].vertex, curve.segm[k].vertex) / denom
			dd = fabs(dd)
			if dd > 1 {
				alpha = 1 - 1.0/dd
			} else {
				alpha = 0
			}
			alpha = alpha / 0.75
		} else {
			alpha = 4 / 3.0
		}
		curve.segm[j].alpha0 = alpha // remember "original" value of alpha

		if alpha >= alphamax { // pointed corner
			curve.segm[j].Type = TypeCorner
			curve.segm[j].Pnt[1] = curve.segm[j].vertex
			curve.segm[j].Pnt[2] = p4
		} else {
			if alpha < 0.55 {
				alpha = 0.55
			} else if alpha > 1 {
				alpha = 1
			}
			p2 := interval(.5+.5*alpha, curve.segm[i].vertex, curve.segm[j].vertex)
			p3 := interval(.5+.5*alpha, curve.segm[k].vertex, curve.segm[j].vertex)
			curve.segm[j].Type = TypeBezier
			curve.segm[j].Pnt[0] = p2
			curve.segm[j].Pnt[1] = p3
			curve.segm[j].Pnt[2] = p4
		}
		curve.segm[j].alpha = alpha // store the "cropped" value of alpha
		curve.segm[j].beta = 0.5
	}
	curve.alphacurve = true
}

// Stage 5: Curve optimization (Sec. 2.4)

type opti struct {
	pen   float64
	c     [2]Point
	t, s  float64
	alpha float64
}

// optiPenalty calculates best fit from i+.5 to j+.5.  Assume i<j (cyclically).
// Return 0 and set badness and parameters (alpha, beta), if possible. Return 1 if impossible.
func optiPenalty(pp *privPath, i, j int, opttolerance float64, convc []int, areac []float64) (res opti, err error) {
	m := len(pp.Curve.segm)

	// check convexity, corner-freeness, and maximum bend < 179 degrees

	if i == j { // sanity - a full loop can never be an opticurve
		err = fmt.Errorf("full loop can never be an opticurve")
		return
	}

	k := i
	i1 := mod(i+1, m)
	k1 := mod(k+1, m)
	conv := convc[k1]
	if conv == 0 {
		err = fmt.Errorf("conv is zero")
		return
	}
	d := ddist(pp.Curve.segm[i].vertex, pp.Curve.segm[i1].vertex)
	for k := k1; k != j; k = k1 {
		k1 = mod(k+1, m)
		k2 := mod(k+2, m)
		if convc[k1] != conv {
			err = fmt.Errorf("convc != conv")
			return
		}
		if int(signf(cprod(pp.Curve.segm[i].vertex, pp.Curve.segm[i1].vertex, pp.Curve.segm[k1].vertex, pp.Curve.segm[k2].vertex))) != conv {
			err = fmt.Errorf("sign 1")
			return
		}
		if iprod1(pp.Curve.segm[i].vertex, pp.Curve.segm[i1].vertex, pp.Curve.segm[k1].vertex, pp.Curve.segm[k2].vertex) < d*ddist(pp.Curve.segm[k1].vertex, pp.Curve.segm[k2].vertex)*cCOS179 {
			err = fmt.Errorf("iprod1")
			return
		}
	}

	// the curve we're working in:
	p0 := pp.Curve.segm[mod(i, m)].Pnt[2]
	p1 := pp.Curve.segm[mod(i+1, m)].vertex
	p2 := pp.Curve.segm[mod(j, m)].vertex
	p3 := pp.Curve.segm[mod(j, m)].Pnt[2]

	// determine its area
	area := areac[j] - areac[i]
	area -= dpara(pp.Curve.segm[0].vertex, pp.Curve.segm[i].Pnt[2], pp.Curve.segm[j].Pnt[2]) / 2
	if i >= j {
		area += areac[m]
	}

	//  	 find intersection o of p0p1 and p2p3. Let t,s such that o =
	// 	 interval(t,p0,p1) = interval(s,p3,p2). Let A be the area of the
	// 	 triangle (p0,o,p3).

	A1 := dpara(p0, p1, p2)
	A2 := dpara(p0, p1, p3)
	A3 := dpara(p0, p2, p3)
	// A4 = dpara(p1, p2, p3)
	A4 := A1 + A3 - A2

	if A2 == A1 { // this should never happen
		err = fmt.Errorf("A2 vs A1")
		return
	}

	t := A3 / (A3 - A4)
	s := A2 / (A2 - A1)
	A := A2 * t / 2.0

	if A == 0.0 { // this should never happen
		err = fmt.Errorf("A == 0")
		return
	}

	R := area / A                   // relative area
	alpha := 2 - math.Sqrt(4-R/0.3) // overall alpha for p0-o-p3 curve

	res.c[0] = interval(t*alpha, p0, p1)
	res.c[1] = interval(s*alpha, p3, p2)
	res.alpha = alpha
	res.t = t
	res.s = s

	p1 = res.c[0]
	p2 = res.c[1] // the proposed curve is now (p0,p1,p2,p3)

	res.pen = 0

	// calculate penalty
	// check tangency with edges
	for k = mod(i+1, m); k != j; k = k1 {
		k1 = mod(k+1, m)
		t = tangent(p0, p1, p2, p3, pp.Curve.segm[k].vertex, pp.Curve.segm[k1].vertex)
		if t < -0.5 {
			err = fmt.Errorf("t < -0.5")
			return
		}
		pt := bezier(t, p0, p1, p2, p3)
		d = ddist(pp.Curve.segm[k].vertex, pp.Curve.segm[k1].vertex)
		if d == 0.0 { // this should never happen
			err = fmt.Errorf("d == 0")
			return
		}
		d1 := dpara(pp.Curve.segm[k].vertex, pp.Curve.segm[k1].vertex, pt) / d
		if fabs(d1) > opttolerance {
			err = fmt.Errorf("tollerance")
			return
		}
		if iprod(pp.Curve.segm[k].vertex, pp.Curve.segm[k1].vertex, pt) < 0 || iprod(pp.Curve.segm[k1].vertex, pp.Curve.segm[k].vertex, pt) < 0 {
			err = fmt.Errorf("prods")
			return
		}
		res.pen += d1 * d1
	}

	// check corners
	for k = i; k != j; k = k1 {
		k1 = mod(k+1, m)
		t = tangent(p0, p1, p2, p3, pp.Curve.segm[k].Pnt[2], pp.Curve.segm[k1].Pnt[2])
		if t < -0.5 {
			err = fmt.Errorf("t < -0.5 s")
			return
		}
		pt := bezier(t, p0, p1, p2, p3)
		d = ddist(pp.Curve.segm[k].Pnt[2], pp.Curve.segm[k1].Pnt[2])
		if d == 0.0 { // this should never happen
			err = fmt.Errorf("d == 0 s")
			return
		}
		d1 := dpara(pp.Curve.segm[k].Pnt[2], pp.Curve.segm[k1].Pnt[2], pt) / d
		d2 := dpara(pp.Curve.segm[k].Pnt[2], pp.Curve.segm[k1].Pnt[2], pp.Curve.segm[k1].vertex) / d
		d2 *= 0.75 * pp.Curve.segm[k1].alpha
		if d2 < 0 {
			d1 = -d1
			d2 = -d2
		}
		if d1 < d2-opttolerance {
			err = fmt.Errorf("tolerance 2")
			return
		}
		if d1 < d2 {
			res.pen += (d1 - d2) * (d1 - d2)
		}
	}

	return
}

// optimize the path p, replacing sequences of Bezier segments by a
// single segment when possible.
func (pp *privPath) opticurve(opttolerance float64) error {
	m := len(pp.Curve.segm)

	var (
		pt    = make([]int, m+1)
		pen   = make([]float64, m+1)
		leng  = make([]int, m+1)
		opt   = make([]opti, m+1)
		convc = make([]int, m)       // pre-computed convexities
		areac = make([]float64, m+1) // cache for fast area computation
	)

	// pre-calculate convexity: +1 = right turn, -1 = left turn, 0 = corner
	for i := 0; i < m; i++ {
		if pp.Curve.segm[i].Type == TypeBezier {
			convc[i] = int(signf(dpara(pp.Curve.segm[mod(i-1, m)].vertex, pp.Curve.segm[i].vertex, pp.Curve.segm[mod(i+1, m)].vertex)))
		} else {
			convc[i] = 0
		}
	}

	// pre-calculate areas
	area := 0.0
	areac[0] = 0.0
	p0 := pp.Curve.segm[0].vertex
	for i := 0; i < m; i++ {
		i1 := mod(i+1, m)
		if pp.Curve.segm[i1].Type == TypeBezier {
			alpha := pp.Curve.segm[i1].alpha
			area += 0.3 * alpha * (4 - alpha) * dpara(pp.Curve.segm[i].Pnt[2], pp.Curve.segm[i1].vertex, pp.Curve.segm[i1].Pnt[2]) / 2
			area += dpara(p0, pp.Curve.segm[i].Pnt[2], pp.Curve.segm[i1].Pnt[2]) / 2
		}
		areac[i+1] = area
	}

	pt[0] = -1
	pen[0] = 0
	leng[0] = 0

	// Fixme: we always start from a fixed point -- should find the best curve cyclically

	for j := 1; j <= m; j++ {
		// calculate best path from 0 to j
		pt[j] = j - 1
		pen[j] = pen[j-1]
		leng[j] = leng[j-1] + 1

		for i := j - 2; i >= 0; i-- {
			o, err := optiPenalty(pp, i, mod(j, m), opttolerance, convc, areac)
			if err != nil {
				break
				//return fmt.Errorf("opti: %s", err)
			}
			if leng[j] > leng[i]+1 || (leng[j] == leng[i]+1 && pen[j] > pen[i]+o.pen) {
				pt[j] = i
				pen[j] = pen[i] + o.pen
				leng[j] = leng[i] + 1
				opt[j] = o
			}
		}
	}
	om := leng[m]
	pp.OCurve = newPrivCurve(om)
	s := make([]float64, om)
	t := make([]float64, om)

	j := m
	for i := om - 1; i >= 0; i-- {
		if pt[j] == j-1 {
			pp.OCurve.segm[i] = pp.Curve.segm[mod(j, m)]
			s[i], t[i] = 1.0, 1.0
		} else {
			pp.OCurve.segm[i] = segment{
				Segment: Segment{
					Type: TypeBezier,
					Pnt:  [3]Point{opt[j].c[0], opt[j].c[1], pp.Curve.segm[mod(j, m)].Pnt[2]},
				},
				vertex: interval(opt[j].s, pp.Curve.segm[mod(j, m)].Pnt[2], pp.Curve.segm[mod(j, m)].vertex),
				alpha:  opt[j].alpha,
				alpha0: opt[j].alpha,
			}
			s[i] = opt[j].s
			t[i] = opt[j].t
		}
		j = pt[j]
	}

	// calculate beta parameters
	for i := 0; i < om; i++ {
		i1 := mod(i+1, om)
		pp.OCurve.segm[i].beta = s[i] / (s[i] + t[i1])
	}
	pp.OCurve.alphacurve = true
	return nil
}

func processPath(plist []Path, param *Params) error {
	// call downstream function with each path
	for i := range plist {
		if err := processOnePath(&plist[i], param); err != nil {
			return err
		}
	}
	return nil
}

func processOnePath(p *Path, param *Params) error {
	p.priv.calcSums()
	p.priv.calcLon()
	p.priv.bestpolygon()
	p.priv.adjustVertices()
	if p.Sign == -1 { // reverse orientation of negative paths
		reverse(&p.priv.Curve)
	}
	smooth(&p.priv.Curve, param.AlphaMax)
	if param.OptiCurve {
		if err := p.priv.opticurve(param.OptTolerance); err != nil {
			return err
		}
		p.priv.FCurve = &p.priv.OCurve
	} else {
		p.priv.FCurve = &p.priv.Curve
	}
	p.Curve = p.priv.FCurve.ToCurve()
	//p.priv = nil
	return nil
}
