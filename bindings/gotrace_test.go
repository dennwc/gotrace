package gotraceb

import (
	"crypto/sha1"
	"fmt"
	gtr "github.com/dennwc/gotrace"
	"image"
	"image/color"
	_ "image/png"
	"os"
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
	data := imageToBitmap(genTestBitmap(), gtr.ThresholdAlpha)
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

func TestTraceBoth(t *testing.T) {
	img := genTestBitmap()
	bm := gtr.NewBitmapFromImage(img, gtr.ThresholdAlpha)
	{
		paths, err := Trace(bm, nil)
		if err != nil {
			t.Fatal(err)
		} else if len(paths) != 1 {
			t.Fatal("wrong paths count:", len(paths))
		} else if paths[0].Sign != +1 {
			t.Fatal("wrong sign in path:", paths[0].Sign)
		}
		t.Log(paths)
	}
	{
		paths, err := Trace(bm, nil)
		if err != nil {
			t.Fatal(err)
		} else if len(paths) != 1 {
			t.Fatal("wrong paths count:", len(paths))
		} else if paths[0].Sign != +1 {
			t.Fatal("wrong sign in path:", paths[0].Sign)
		}
		t.Log(paths)
	}
}

func loadBitmap(t *testing.T) *gtr.Bitmap {
	file, err := os.Open("../../inp.png")
	if err != nil {
		t.Skip(err)
	}
	defer file.Close()
	img, _, err := image.Decode(file)
	if err != nil {
		t.Fatal(err)
	}
	bm := gtr.NewBitmapFromImage(img, func(x, y int, c color.Color) bool {
		r, _, _, _ := c.RGBA()
		return (float64(r)/0xffff)*0xff > 128
	})
	return bm
}

func TestTraceCompare(t *testing.T) {
	bm := loadBitmap(t)

	hashs := func(o interface{}) string {
		h := sha1.New()
		fmt.Fprint(h, o)
		return fmt.Sprintf("%x", h.Sum(nil))
	}

	paths1, err := gtr.Trace(bm, nil)
	if err != nil {
		t.Fatal(err)
	}

	paths2, err := Trace(bm, nil)
	if err != nil {
		t.Fatal(err)
	}

	hmap := make(map[string]bool)
	for _, p := range paths1 {
		hmap[hashs(p)] = true
	}
	for _, p := range paths2 {
		h := hashs(p)
		_, ok := hmap[h]
		if ok {
			delete(hmap, h)
		} else {
			hmap[h] = false
		}
	}
	if len(hmap) > 0 {
		for _, p := range paths1 {
			if _, ok := hmap[hashs(p)]; ok {
				t.Log(p)
			}
		}
		for _, p := range paths2 {
			if _, ok := hmap[hashs(p)]; ok {
				t.Log(p)
			}
		}
		t.Fatal("different paths", hmap)
	}

	if len(paths1) != len(paths2) {
		t.Fatalf("paths count differs: %d vs %d", len(paths1), len(paths2))
	}
	if hashs(paths1) != hashs(paths2) {
		t.Log("output differs")
	}
}
