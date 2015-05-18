# Interpolations

[![Build Status](https://travis-ci.org/tlycken/Interpolations.jl.svg?branch=master)](https://travis-ci.org/tlycken/Interpolations.jl)

This package is the continuation of Tim Holy's work on [Grid.jl](https://github.com/timholy/Grid.jl), with the explicit goal of reaching feature parity for [B-splines](https://en.wikipedia.org/wiki/B-spline) in the near future (PRs are more than welcome!).

## Installation

Just

```
Pkg.add("Interpolations")
```

from the Julia REPL.

## Usage

*Note: the API might change in future versions, if necessary to accomodate new functionality.*

```
using Interpolations
x = 1:5 # must be a unit range starting at 1
coarse = sin(x)
itp1 = Interpolation(coarse, Linear(OnGrid()), ExtrapPeriodic())
sinx = [itp1[x] for x in -pi:.01:3pi]
cosx = [gradient(itp,x)[1] for x in -pi:.1:3pi]
```

## Options

To interpolate this data using `Interpolations.jl` you provide the constructor with a few arguments that describe *how* the interpolation should be performed.

### Interpolation

The first argument describes the interpolation type in terms of *degree*, *grid behavior* and, if necessary, *boundary conditions*. There are currently three degrees available: `Constant`, `Linear` and `Quadratic`, corresponding to B-splines of degree 0, 1 and 2, respectively.

You also have to specify what *grid representation* you want. There are currently two choices: `OnGrid`, in which the data points are assumed to lie *on* the boundaries of the interpolation interval, and `OnCell` in which the data points are assumed to lie on half-intervals between cell boundaries.

B-splines of quadratic or higher degree require solving an equation system to obtain the interpolation coefficients, and for that you must specify a *boundary condition* that is applied to close the system. The following boundary conditions are implemented: `Flat`, `Line`, `Free`, `Periodic` and `Reflect`; their mathematical implications are described in detail in the pdf document under /doc/latex.

### Extrapolation

The final argument to the `Interpolation` constructor describes what happens if you try to index into the interpolation object with coordinates outside of `[1, size(data,d)]` in any dimension `d`. The implemented boundary conditions are `ExtrapError`, `ExtrapNaN`, `ExtrapConstant`, `ExtrapLinear`, `ExtrapPeriodic` and `ExtrapReflect`.

Note that for extrapolating a periodic data set in a way that makes sense, you should specify *both* the `Periodic` boundary condition *and* the `ExtrapPeriodic` extrapolation behavior, and similarly for reflecting extrapolation. There might be API changes in the future to make these behaviors easier to specify.

## More examples

There's an [IJulia notebook](http://nbviewer.ipython.org/github/tlycken/Interpolations.jl/blob/master/doc/Interpolations.jl.ipynb) that shows off some of the current functionality, and outlines where we're headed. I try to keep it up to date when I make any significant improvements and/or breaking changes, but if it's not, do file a PR.

## Contributing

Work is very much in progress, but and help is always welcome. If you want to help out but don't know where to start, take a look at issue [#5 - our feature wishlist](https://github.com/tlycken/Interpolations.jl/issues/5) =)

Contributions in any form are appreciated, but my time is limited and your code is most likely to be included quickly in the library if it comes in a pull request with at least a few sanity tests (see the `/test` directory for some examples).
