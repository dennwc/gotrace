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
	TypeBezier = SegmType(1)
	TypeCorner = SegmType(2)
)

type SegmType int

const (
	TurnBlack    = TurnPolicy(0)
	TurnWhite    = TurnPolicy(1)
	TurnLeft     = TurnPolicy(2)
	TurnRight    = TurnPolicy(3)
	TurnMinority = TurnPolicy(4)
	TurnMajority = TurnPolicy(5)
	TurnRandom   = TurnPolicy(6)
)

type TurnPolicy int

type Point struct {
	X, Y float64
}

type Segment struct {
	Type SegmType
	Pnt  [3]Point
}

type Path struct {
	Area  int
	Sign  int
	Curve []Segment

	Childs []Path

	priv *privPath
}

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

type Params struct {
	TurdSize     int
	TurnPolicy   TurnPolicy
	AlphaMax     float64
	OptiCurve    bool
	OptTolerance float64
}

var Defaults = Params{
	TurdSize:     2,
	TurnPolicy:   TurnMinority,
	AlphaMax:     1.0,
	OptiCurve:    true,
	OptTolerance: 0.2,
}

func Trace(bm *Bitmap, param *Params) ([]Path, error) {
	if param == nil {
		param = &Defaults
	}
	return bm.toPathList(param)
}

func ThresholdAlpha(x, y int, cl color.Color) bool {
	_, _, _, a := cl.RGBA()
	return a != 0
}

func NewBitmapFromImage(img image.Image, threshold func(x, y int, cl color.Color) bool) *Bitmap {
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
