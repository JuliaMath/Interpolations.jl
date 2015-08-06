# Interpolations

[![Build Status](https://travis-ci.org/tlycken/Interpolations.jl.svg?branch=master)](https://travis-ci.org/tlycken/Interpolations.jl)

This package implements a variety of interpolation schemes for the
Julia langauge.  It has the goals of ease-of-use, broad algorithmic
support, and exceptional performance.

This package is still relatively new. Currently its support is best
for [B-splines](https://en.wikipedia.org/wiki/B-spline), but the API
has been designed with intent to support more options. Pull-requests
are more than welcome!  It should be noted that the API may continue
to evolve over time.

Other interpolation packages for Julia include:
- [Grid.jl](https://github.com/timholy/Grid.jl) (the predecessor of this package)
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
- [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
- [ApproXD.jl](https://github.com/floswald/ApproXD.jl)

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

## Installation

Just

```
Pkg.add("Interpolations")
```

from the Julia REPL.

## General usage

Given an `AbstractArray` `A`, construct an "interpolation object" `itp` as
```jl
itp = interpolate(A, options...)
```
where `options...` (discussed below) controls the type of
interpolation you want to perform.  To evaluate the interpolation at
position `(x, y, ...)`, simply do
```jl
v = itp[x, y, ...]
```

Some interpolation objects support computation of the gradient, which
can be obtained as
```jl
g = gradient(itp, x, y, ...)
```
or, if you're evaluating the gradient repeatedly, a somewhat more
efficient option is
```jl
gradient!(g, itp, x, y, ...)
```
where `g` is a pre-allocated vector.

`A` may have any element type that supports the operations of addition
and multiplication.  Examples include scalars like `Float64`, `Int`,
and `Rational`, but also multi-valued types like `RGB` color vectors.

Positions `(x, y, ...)` are n-tuples of numbers. Typically these will
be real-valued (not necessarily integer-valued), but can also be of types
such as [DualNumbers](https://github.com/JuliaDiff/DualNumbers.jl) if
you want to verify the computed value of gradients.


## Control of interpolation algorithm

### B-splines

The interpolation type is described in terms of *degree*, *grid behavior* and, if necessary, *boundary conditions*. There are currently three degrees available: `Constant`, `Linear` and `Quadratic`, corresponding to B-splines of degree 0, 1 and 2, respectively.

You also have to specify what *grid representation* you want. There are currently two choices: `OnGrid`, in which the supplied data points are assumed to lie *on* the boundaries of the interpolation interval, and `OnCell` in which the data points are assumed to lie on half-intervals between cell boundaries.

B-splines of quadratic or higher degree require solving an equation system to obtain the interpolation coefficients, and for that you must specify a *boundary condition* that is applied to close the system. The following boundary conditions are implemented: `Flat`, `Line` (alternatively, `Natural`), `Free`, `Periodic` and `Reflect`; their mathematical implications are described in detail in the pdf document under `/doc/latex`.

Some examples:
```jl
# Nearest-neighbor interpolation
itp = interpolate(a, BSpline{Constant}, OnCell)
v = itp[5.4]   # returns a[5]

# (Multi)linear interpolation
itp = interpolate(A, BSpline{Linear}, OnGrid)
v = itp[3.2, 4.1]  # returns 0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5])

# Quadratic interpolation with reflecting boundary conditions
# Quadratic is the lowest order that has continuous gradient
itp = interpolate(A, BSpline{Quadratic{Reflect}}, OnCell)

# Linear interpolation in the first dimension, and no interpolation (just lookup) in the second
itp = interpolate(A, Tuple{BSpline{Linear}, BSpline{NoInterp}}, OnGrid)
v = itp[3.65, 5]  # returns  0.35*A[3,5] + 0.65*A[4,5]
```
There are more options available, for example:
```
# In-place interpolation
itp = interpolate!(A, BSpline{Quadratic{InPlace}}, OnCell)
```
which destroys the input `A` but also does not need to allocate as much memory.

## Extrapolation

The call to `extrapolate` defines what happens if you try to index into the interpolation object with coordinates outside of `[1, size(data,d)]` in any dimension `d`. The implemented boundary conditions are `Throw` and `Flat`, with more options planned.

## More examples

There's an [IJulia notebook](http://nbviewer.ipython.org/github/tlycken/Interpolations.jl/blob/master/doc/Interpolations.jl.ipynb) that shows off some of the current functionality, and outlines where we're headed. I try to keep it up to date when I make any significant improvements and/or breaking changes, but if it's not, do file a PR.

## Contributing

Work is very much in progress, but and help is always welcome. If you want to help out but don't know where to start, take a look at issue [#5 - our feature wishlist](https://github.com/tlycken/Interpolations.jl/issues/5) =) There is also some [developer documentation](doc/devdocs.md) that may help you understand how things work internally.

Contributions in any form are appreciated, but the best pull requests come with tests!
