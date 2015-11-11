package gotrace

import (
	"image"
	"image/color"
	"math/rand"
	"testing"
)

func TestBuild(t *testing.T) {}

func genTestBitmap() image.Image {
	img := image.NewRGBA(image.Rect(0, 0, 12, 12))
	for y := 0; y < 12; y++ {
		for x := 0; x < 12; x++ {
			if x >= 11-y {
				img.Set(x, y, color.Black)
			} else {
				img.Set(x, y, color.Transparent)
			}
		}
	}
	return img
}

func TestBitmap(t *testing.T) {
	bm := NewBitmap(1543, 1234)
	i := 0
	check := make([]bool, bm.W*bm.H)
	for y := 0; y < bm.H; y++ {
		for x := 0; x < bm.W; x++ {
			v := rand.Intn(10) > 4
			check[i] = v
			bm.Set(x, y, v)
			i++
		}
	}
	i = 0
	for y := 0; y < bm.H; y++ {
		for x := 0; x < bm.W; x++ {
			if check[i] != bm.Get(x, y) {
				t.Fatal("failed")
			}
			i++
		}
	}
}

func TestTrace(t *testing.T) {
	bm := NewBitmapFromImage(genTestBitmap(), ThresholdAlpha)
	paths, err := Trace(bm, nil)
	if err != nil {
		t.Fatal(err)
	} else if len(paths) != 1 {
		t.Fatal("wrong paths count:", len(paths))
	} else if paths[0].Sign != +1 {
		t.Fatal("wrong sign in path:", paths[0].Sign)
	} else if paths[0].Area != 78 {
		t.Fatal("wrong area of path:", paths[0].Area)
	}
}

func BenchmarkBitmap(b *testing.B) {
	bm := NewBitmap(1543, 1234)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x, y, v := i%bm.W, i%bm.H, i%2 == 0
		for j := 0; j < 1000; j++ {
			bm.Set(x, y, v)
			bm.Get(x, y)
		}
	}
}
