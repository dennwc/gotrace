// Package gotrace provides bindings for potracelib - bitmap to vector graphics converter.
// More info at http://potrace.sourceforge.net/potracelib.pdf
package gotrace

// #cgo LDFLAGS: -lpotrace
// #include <stdlib.h>
// #include <potracelib.h>
// int wordSize() { return sizeof(potrace_word); }
// potrace_param_t* newParams() { return malloc(sizeof(potrace_param_t)); }
// potrace_bitmap_t* newBitmap() { return malloc(sizeof(potrace_bitmap_t)); }
import "C"

import (
	"bytes"
	"fmt"
	"image"
	"image/color"
	"io"
	"unsafe"
)

var wordSize = int(C.wordSize()) * 8

type TurnPolicy int

const (
	TurnPolicyBlack    = TurnPolicy(C.POTRACE_TURNPOLICY_BLACK)
	TurnPolicyWhite    = TurnPolicy(C.POTRACE_TURNPOLICY_WHITE)
	TurnPolicyLeft     = TurnPolicy(C.POTRACE_TURNPOLICY_LEFT)
	TurnPolicyRight    = TurnPolicy(C.POTRACE_TURNPOLICY_RIGHT)
	TurnPolicyMinority = TurnPolicy(C.POTRACE_TURNPOLICY_MINORITY)
	TurnPolicyMajority = TurnPolicy(C.POTRACE_TURNPOLICY_MAJORITY)
	TurnPolicyRandom   = TurnPolicy(C.POTRACE_TURNPOLICY_RANDOM)
)

type SegType int

const (
	Bezier = SegType(C.POTRACE_CURVETO)
	Corner = SegType(C.POTRACE_CORNER)
)

func Version() string {
	return C.GoString(C.potrace_version())
}

type Path struct {
	Area  int // enclosed area
	Sign  int // +1 or -1
	Curve []Segment
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
			fmt.Fprintf(buf, "M%f,%f ", c[l].Pt[2][0], c[l].Pt[2][1]) // end point == start point
		}
		if p.Type == Corner {
			fmt.Fprintf(buf, "L%f,%f", p.Pt[1][0], p.Pt[1][1]) // no last point for now - may be closed with Z
		} else if p.Type == Bezier {
			fmt.Fprintf(buf, "C%f,%f %f,%f %f,%f", p.Pt[0][0], p.Pt[0][1], p.Pt[1][0], p.Pt[1][1], p.Pt[2][0], p.Pt[2][1])
		}
		if i == l {
			fmt.Fprintf(buf, " Z")
		} else if p.Type == Corner {
			fmt.Fprintf(buf, " %f,%f ", p.Pt[2][0], p.Pt[2][1]) // last point for corner
		}
	}
	return buf.String()
}

type Segment struct {
	Type SegType
	Pt   [3]Point
}

type Point [2]float64

type Threshold func(color.Color) bool

// Params is a structure to hold tracing parameters
type Params struct {
	TurdSize     int        // area of largest path to be ignored
	TurnPolicy   TurnPolicy // resolves ambiguous turns in path decomposition
	AlphaMax     float64    // corner threshold
	OptCurve     bool       // use curve optimization?
	OptTolerance float64    // curve optimization tolerance

	ThresholdFunc Threshold // function to convert pixels to bits
}

func ThresholdAplha(c color.Color) bool {
	_, _, _, a := c.RGBA()
	return a != 0
}

func Defaults() *Params {
	out := new(Params)
	p := C.potrace_param_default()
	defer C.free(unsafe.Pointer(p))
	out.TurdSize = int(p.turdsize)
	out.TurnPolicy = TurnPolicy(p.turnpolicy)
	out.AlphaMax = float64(p.alphamax)
	if int(p.opticurve) != 0 {
		out.OptCurve = true
	}
	out.OptTolerance = float64(p.opttolerance)

	out.ThresholdFunc = ThresholdAplha
	return out
}

func imageToBitmap(img image.Image, th Threshold) []uint {
	w, h := img.Bounds().Dx(), img.Bounds().Dy()
	ww := w / wordSize
	if w%wordSize != 0 {
		ww++
	}
	words := make([]uint, ww*h)
	for y := 0; y < img.Bounds().Dy(); y++ {
		for x := 0; x < img.Bounds().Dx(); x++ {
			if th(img.At(x, y)) {
				i := y*ww + x/wordSize
				words[i] = words[i] | (1 << uint(wordSize-1-x%wordSize))
			}
		}
	}
	return words
}

func Trace(img image.Image, params *Params) ([]Path, error) {
	if params == nil {
		params = Defaults()
	}
	data := imageToBitmap(img, params.ThresholdFunc)
	w, h := img.Bounds().Dx(), img.Bounds().Dy()
	dy := w / wordSize
	if w%wordSize != 0 {
		dy++
	}
	par := C.newParams()
	defer C.free(unsafe.Pointer(par))
	par.turdsize = C.int(params.TurdSize)
	par.turnpolicy = C.int(params.TurnPolicy)
	par.alphamax = C.double(params.AlphaMax)
	if params.OptCurve {
		par.opticurve = C.int(1)
	} else {
		par.opticurve = C.int(0)
	}
	par.opttolerance = C.double(params.OptTolerance)

	bmp := C.newBitmap()
	defer C.free(unsafe.Pointer(bmp))
	bmp._map = (*C.potrace_word)(unsafe.Pointer(&data[0]))
	bmp.w = C.int(w)
	bmp.h = C.int(h)
	bmp.dy = C.int(dy)

	st, err := C.potrace_trace(par, bmp)
	defer C.potrace_state_free(st)
	if err != nil {
		return nil, err
	} else if int(st.status) != 0 {
		return nil, fmt.Errorf("trace failed")
	}
	var out []Path
	cur := st.plist
	for cur != nil {
		var p Path
		p.Area = int(cur.area)
		if byte(cur.sign) == byte('+') {
			p.Sign = +1
		} else if byte(cur.sign) == byte('-') {
			p.Sign = -1
		}
		// process curve
		p.Curve = make([]Segment, int(cur.curve.n))
		for i := range p.Curve {
			p.Curve[i] = Segment{
				Type: SegType(*(*C.int)(unsafe.Pointer(uintptr(unsafe.Pointer(cur.curve.tag)) + unsafe.Sizeof(C.int(0))*uintptr(i)))),
				Pt: [3]Point{
					pointFromArr(cur.curve.c, i, 0),
					pointFromArr(cur.curve.c, i, 1),
					pointFromArr(cur.curve.c, i, 2),
				},
			}
		}
		out = append(out, p)
		cur = cur.next
	}
	return out, nil
}

func pointFromArr(c *[3]C.potrace_dpoint_t, i, j int) Point {
	p := *(*[3]C.potrace_dpoint_t)(unsafe.Pointer(uintptr(unsafe.Pointer(c)) + unsafe.Sizeof([3]C.potrace_dpoint_t{})*uintptr(i)))
	return Point{float64(p[j].x), float64(p[j].y)}
}

func WriteSvg(w io.Writer, rect image.Rectangle, paths []Path) error {
	fmt.Fprintf(w, `<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20010904//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg version="1.0" xmlns="http://www.w3.org/2000/svg" width="%dpt" height="%dpt" viewBox="0 0 %d %d" preserveAspectRatio="xMidYMid meet">
<g transform="translate(0.0,%d.0) scale(1.0,-1.0)" fill="#000000" stroke="none">%s`,
		rect.Dx(), rect.Dy(), rect.Dx(), rect.Dy(), rect.Dy(), "\n")
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
