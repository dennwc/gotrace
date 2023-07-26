# GoTrace

**WARNING:** This library is **deprecated** in favor of [gotranspile/gotrace](https://github.com/gotranspile/gotrace). It will not receive any updates.

Pure Go implementation of [Potrace](http://potrace.sourceforge.net/potracelib.pdf) vectorization library.
Supports simple SVG output generation.

Also includes [cgo bindings](./bindings) for the original Potrace library.

**Original image**

![Original](http://potrace.sourceforge.net/img/stanford-orig2.png)

**Vectorized image**

![Vectorized](http://potrace.sourceforge.net/img/stanford-smooth2.png)

# Installation
```
go get github.com/dennwc/gotrace
```

# Usage

Process image with an alpha channel, generate SVG:

```Go
bm := gotrace.NewBitmapFromImage(img, nil)
paths, _ := gotrace.Trace(bm, nil)
gotrace.WriteSvg(file, img.Bounds(), paths, "")
```

Custom threshold function:

```Go
bm := gotrace.NewBitmapFromImage(img, func(x, y int, c color.Color) bool {
    r, g, b, _ := c.RGBA()
    return r + g + b > 128
})
paths, _ := gotrace.Trace(bm, nil)
```
