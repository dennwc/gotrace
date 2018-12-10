// Package gotrace is a pure Go implementation of potrace vectorization library.
// More info at http://potrace.sourceforge.net/potracelib.pdf
package gotrace

import (
	"bytes"
	"fmt"
	"image"
	"image/color"
	"io"
)

const (
	// TypeBezier represents Bezier curve segment.
	TypeBezier = SegmType(1)
	// TypeCorner represents a corner segment.
	TypeCorner = SegmType(2)
)

// SegmType is a type of a path segment.
type SegmType int

// Supported turn policies.
const (
	TurnBlack    = TurnPolicy(0)
	TurnWhite    = TurnPolicy(1)
	TurnLeft     = TurnPolicy(2)
	TurnRight    = TurnPolicy(3)
	TurnMinority = TurnPolicy(4)
	TurnMajority = TurnPolicy(5)
	TurnRandom   = TurnPolicy(6)
)

// TurnPolicy is a policy for tracing the bitmap.
type TurnPolicy int

// Point is a point on a 2D plane.
type Point struct {
	X, Y float64
}

// Segment is a pth segment.
type Segment struct {
	Type SegmType
	Pnt  [3]Point
}

// Path represents a curve path.
type Path struct {
	Area  int
	Sign  int
	Curve []Segment

	Childs []Path

	priv *privPath
}

// ToSvgPath converts the path to an SVG path string.
func (p Path) ToSvgPath() string {
	c := p.Curve
	if len(c) == 0 {
		return ""
	}
	l := len(c) - 1
	buf := bytes.NewBuffer(nil)
	for i, p := range c {
		if i == 0 {
			fmt.Fprintf(buf, "M%f,%f ", c[l].Pnt[2].X, c[l].Pnt[2].Y) // end point == start point
		}
		if p.Type == TypeCorner {
			fmt.Fprintf(buf, "L%f,%f", p.Pnt[1].X, p.Pnt[1].Y) // no last point for now - may be closed with Z
		} else if p.Type == TypeBezier {
			fmt.Fprintf(buf, "C%f,%f %f,%f %f,%f", p.Pnt[0].X, p.Pnt[0].Y, p.Pnt[1].X, p.Pnt[1].Y, p.Pnt[2].X, p.Pnt[2].Y)
		}
		if i == l {
			fmt.Fprintf(buf, " Z")
		} else if p.Type == TypeCorner {
			fmt.Fprintf(buf, " %f,%f ", p.Pnt[2].X, p.Pnt[2].Y) // last point for corner
		}
	}
	return buf.String()
}

// Params is a set of tracing parameters.
type Params struct {
	TurdSize     int
	TurnPolicy   TurnPolicy
	AlphaMax     float64
	OptiCurve    bool
	OptTolerance float64
}

// Defaults cores default parameters for tracing.
var Defaults = Params{
	TurdSize:     2,
	TurnPolicy:   TurnMinority,
	AlphaMax:     1.0,
	OptiCurve:    true,
	OptTolerance: 0.2,
}

// Trace the bitmap using specified parameters for the algorithm.
// If parameters is nil, defaults will be used.
func Trace(bm *Bitmap, param *Params) ([]Path, error) {
	if param == nil {
		param = &Defaults
	}
	return bm.toPathList(param)
}

// ThresholdFunc is a function that converts a color to a a single bit value.
type ThresholdFunc func(x, y int, cl color.Color) bool

// ThresholdAlpha is a threshold function that cuts fully transparent parts from the image.
func ThresholdAlpha(_, _ int, cl color.Color) bool {
	_, _, _, a := cl.RGBA()
	return a != 0
}

// NewBitmapFromImage converts an image to a bitmap using a threshold function.
// If function is not specified, an Alpha channel will be used for the threshold.
func NewBitmapFromImage(img image.Image, threshold ThresholdFunc) *Bitmap {
	if threshold == nil {
		threshold = ThresholdAlpha
	}
	rect := img.Bounds()
	bm := NewBitmap(rect.Dx(), rect.Dy())
	for y := 0; y < rect.Dy(); y++ {
		for x := 0; x < rect.Dx(); x++ {
			bm.Set(x, y, threshold(x, y, img.At(x, y)))
		}
	}
	return bm
}

// WriteSvg writes an SVG file with specified paths. Dimensions of an image will be set to rect,
// and paths will be of a specified color (#rrggbb).
func WriteSvg(w io.Writer, rect image.Rectangle, paths []Path, color string) error {
	if color == "" {
		color = "#000000"
	}
	fmt.Fprintf(w, `<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20010904//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg version="1.0" xmlns="http://www.w3.org/2000/svg" width="%dpt" height="%dpt" viewBox="0 0 %d %d" preserveAspectRatio="xMidYMid meet">
<g fill="%s" stroke="none">%s`,
		rect.Dx(), rect.Dy(), rect.Dx(), rect.Dy(), color, "\n")
	for i, p := range paths {
		if p.Sign > 0 {
			fmt.Fprint(w, `<path d="`)
		}
		fmt.Fprint(w, p.ToSvgPath())
		if i+1 == len(paths) || paths[i+1].Sign > 0 {
			fmt.Fprintln(w, `"/>`)
		} else {
			fmt.Fprint(w, " ")
		}
	}
	_, err := fmt.Fprintln(w, `</g></svg>`)
	return err
}
