// Package gotrace provides bindings for potracelib - bitmap to vector graphics converter.
// More info at http://potrace.sourceforge.net/potracelib.pdf
package gotraceb

// #cgo LDFLAGS: -lm
// #include <stdlib.h>
// #include "bitmap.h"
// #include "potracelib.h"
// int wordSize() { return sizeof(potrace_word); }
// potrace_param_t* newParams() { return malloc(sizeof(potrace_param_t)); }
// static inline void bmSet(potrace_bitmap_t *bm, int x, int y) { BM_SET(bm,x,y); }
// static inline void bmClr(potrace_bitmap_t *bm, int x, int y) { BM_CLR(bm,x,y); }
// static inline int bmGet(potrace_bitmap_t *bm, int x, int y) { return BM_GET(bm,x,y); }
// static inline potrace_word bmRaw(potrace_bitmap_t *bm, int i) { return bm->map[i]; }
import "C"

import (
	//"bytes"
	"fmt"
	"image"
	//"io"
	gtr "github.com/dennwc/gotrace"
	"image/color"
	"unsafe"
)

var wordSize = int(C.wordSize()) * 8

func Version() string {
	return C.GoString(C.potrace_version())
}

/*func (p Path) ToSvgPath() string {
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
}*/

func imageToBitmap(img image.Image, th func(x, y int, c color.Color) bool) []uint {
	w, h := img.Bounds().Dx(), img.Bounds().Dy()
	ww := w / wordSize
	if w%wordSize != 0 {
		ww++
	}
	words := make([]uint, ww*h)
	for y := 0; y < img.Bounds().Dy(); y++ {
		for x := 0; x < img.Bounds().Dx(); x++ {
			if th(x, y, img.At(x, y)) {
				i := y*ww + x/wordSize
				words[i] = words[i] | (1 << uint(wordSize-1-x%wordSize))
			}
		}
	}
	return words
}

func Trace(bm *gtr.Bitmap, param *gtr.Params) ([]gtr.Path, error) {
	if param == nil {
		param = &gtr.Defaults
	}
	bmc := C.bm_new(C.int(bm.W), C.int(bm.H))
	defer C.bm_free(bmc)
	for y := 0; y < bm.H; y++ {
		for x := 0; x < bm.W; x++ {
			if bm.Get(x, y) {
				C.bmSet(bmc, C.int(x), C.int(y))
			} else {
				C.bmClr(bmc, C.int(x), C.int(y))
			}
		}
	}

	par := C.newParams()
	defer C.free(unsafe.Pointer(par))
	par.turdsize = C.int(param.TurdSize)
	par.turnpolicy = C.int(param.TurnPolicy)
	par.alphamax = C.double(param.AlphaMax)
	if param.OptiCurve {
		par.opticurve = C.int(1)
	} else {
		par.opticurve = C.int(0)
	}
	par.opttolerance = C.double(param.OptTolerance)

	st, err := C.potrace_trace(par, bmc)
	defer C.potrace_state_free(st)
	if err != nil {
		return nil, err
	} else if int(st.status) != 0 {
		return nil, fmt.Errorf("trace failed")
	}
	return convPaths(st.plist), nil
}

func convPath(cur *C.potrace_path_t) (p gtr.Path) {
	p.Area = int(cur.area)
	if byte(cur.sign) == byte('+') {
		p.Sign = +1
	} else if byte(cur.sign) == byte('-') {
		p.Sign = -1
	}
	// process curve
	p.Curve = make([]gtr.Segment, int(cur.curve.n))
	for i := range p.Curve {
		p.Curve[i] = gtr.Segment{
			Type: gtr.SegmType(*(*C.int)(unsafe.Pointer(uintptr(unsafe.Pointer(cur.curve.tag)) + unsafe.Sizeof(C.int(0))*uintptr(i)))),
			Pnt: [3]gtr.Point{
				pointFromArr(cur.curve.c, i, 0),
				pointFromArr(cur.curve.c, i, 1),
				pointFromArr(cur.curve.c, i, 2),
			},
		}
	}
	return
}

func convPaths(cur *C.potrace_path_t) (out []gtr.Path) {
	if cur == nil {
		return nil
	}
	for ; cur != nil; cur = cur.sibling {
		p := convPath(cur)
		p.Childs = convPaths(cur.childlist)
		out = append(out, p)
	}
	return
}

func pointFromArr(c *[3]C.potrace_dpoint_t, i, j int) gtr.Point {
	p := *(*[3]C.potrace_dpoint_t)(unsafe.Pointer(uintptr(unsafe.Pointer(c)) + unsafe.Sizeof([3]C.potrace_dpoint_t{})*uintptr(i)))
	return gtr.Point{float64(p[j].x), float64(p[j].y)}
}

/*func WriteSvg(w io.Writer, rect image.Rectangle, paths []Path, color string) error {
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
}*/
