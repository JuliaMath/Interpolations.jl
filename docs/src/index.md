
# Interpolations

[![Build Status](https://travis-ci.org/JuliaMath/Interpolations.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Interpolations.jl)
[![PkgEval Status](http://pkg.julialang.org/badges/Interpolations_0.4.svg)](http://pkg.julialang.org/?pkg=Interpolations)
[![Interpolations](http://pkg.julialang.org/badges/Interpolations_0.5.svg)](http://pkg.julialang.org/?pkg=Interpolations)

**NEWS** v0.9 was a breaking release. See the [news](../../NEWS.md) for details on how to update.

This package implements a variety of interpolation schemes for the
Julia language.  It has the goals of ease-of-use, broad algorithmic
support, and exceptional performance.

Currently this package's support is best
for [B-splines](https://en.wikipedia.org/wiki/B-spline) and also
supports irregular grids.  However, the API has been designed with
intent to support more options. Pull-requests are more than welcome!
It should be noted that the API may continue to evolve over time.

Other interpolation packages for Julia include:
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
- [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
- [ApproXD.jl](https://github.com/floswald/ApproXD.jl)

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

## Installation

Just

```
using Pkg
Pkg.add("Interpolations")
```

from the Julia REPL.
