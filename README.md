# GoTrace
Go bindings for [Potrace](http://potrace.sourceforge.net/potracelib.pdf) vectorization library.
Supports simple SVG output generation.

**Original image**

![Original](http://potrace.sourceforge.net/img/stanford-orig2.png)

**Vectorized image**

![Vectorized](http://potrace.sourceforge.net/img/stanford-smooth2.png)

# Installation
```
apt-get install libpotrace-dev
go get github.com/dennwc/gotrace
```

# Usage
Process image, generate SVG:
```Go
paths, _ := gotrace.Trace(img, nil)
gotrace.WriteSvg(file, img.Bounds(), paths)
```

Custom threshold function:
```Go
params := gotrace.Defaults()
params.ThresholdFunc = func(c color.Color) bool {
  r,g,b,_ := c.RGBA()
  return r+g+b > 128
}
paths, _ := gotrace.Trace(img, params)
```
