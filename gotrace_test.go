package gotrace

import (
	"image"
	"image/color"
	"testing"
	"unsafe"
)

func TestWordSize(t *testing.T) {
	if uintptr(wordSize) != unsafe.Sizeof(uint(0))*8 {
		t.Fatalf("%d vs %s", wordSize, unsafe.Sizeof(uint(0))*8)
	}
}

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

func TestPotraceBitmap(t *testing.T) {
	data := imageToBitmap(genTestBitmap(), ThresholdAplha)
	if len(data) != 12 {
		t.Fatal("wrong size:", len(data))
	}
	start := uint(0x00100000)
	if wordSize == 64 {
		start = uint(0x0010000000000000)
	}
	last := start
	for i := range data {
		if data[i] != last {
			t.Fatalf("wrong representation (%d): 0x%016x vs 0x%016x", i, data[i], last)
		}
		last = (last << 1) | start
	}
}

func TestTrace(t *testing.T) {
	img := genTestBitmap()
	paths, err := Trace(img, nil)
	if err != nil {
		t.Fatal(err)
	} else if len(paths) != 1 {
		t.Fatal("wrong paths count:", len(paths))
	} else if paths[0].Sign != +1 {
		t.Fatal("wrong sign in path:", paths[0].Sign)
	}
}
