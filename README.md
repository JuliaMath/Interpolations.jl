# Interpolations

[![Build Status](https://travis-ci.org/JuliaMath/Interpolations.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Interpolations.jl)
[![PkgEval Status](http://pkg.julialang.org/badges/Interpolations_0.4.svg)](http://pkg.julialang.org/?pkg=Interpolations)
[![Interpolations](http://pkg.julialang.org/badges/Interpolations_0.5.svg)](http://pkg.julialang.org/?pkg=Interpolations)

This package implements a variety of interpolation schemes for the
Julia langauge.  It has the goals of ease-of-use, broad algorithmic
support, and exceptional performance.

This package is still relatively new. Currently its support is best
for [B-splines](https://en.wikipedia.org/wiki/B-spline) and also
supports irregular grids.  However, the API has been designed with
intent to support more options. Pull-requests are more than welcome!
It should be noted that the API may continue to evolve over time.

Other interpolation packages for Julia include:
- [Grid.jl](https://github.com/timholy/Grid.jl) (the predecessor of this package)
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
- [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
- [ApproXD.jl](https://github.com/floswald/ApproXD.jl)

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

At the bottom of this page, you can find a "performance shootout"
among these methods (as well as SciPy's `RegularGridInterpolator`).

## Installation

Just

```
Pkg.add("Interpolations")
```

from the Julia REPL.

## General usage

Note: the current version of `Interpolations` supports interpolation evaluation using index calls `[]`, but this feature will be deprecated in future. We highly recommend function calls with `()` as follows.

Given an `AbstractArray` `A`, construct an "interpolation object" `itp` as
```jl
itp = interpolate(A, options...)
```
where `options...` (discussed below) controls the type of
interpolation you want to perform.  This syntax assumes that the
samples in `A` are equally-spaced.

To evaluate the interpolation at position `(x, y, ...)`, simply do
```jl
v = itp(x, y, ...)
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

Some interpolation objects support computation of the hessian, which
can be obtained as
```jl
h = hessian(itp, x, y, ...)
```
or, if you're evaluating the hessian repeatedly, a somewhat more
efficient option is
```jl
hessian!(h, itp, x, y, ...)
```
where `h` is a pre-allocated matrix.

`A` may have any element type that supports the operations of addition
and multiplication.  Examples include scalars like `Float64`, `Int`,
and `Rational`, but also multi-valued types like `RGB` color vectors.

Positions `(x, y, ...)` are n-tuples of numbers. Typically these will
be real-valued (not necessarily integer-valued), but can also be of types
such as [DualNumbers](https://github.com/JuliaDiff/DualNumbers.jl) if
you want to verify the computed value of gradients. You can also use
Julia's iterator objects, e.g.,

```jl
function ongrid!(dest, itp)
    for I in CartesianRange(size(itp))
        dest[I] = itp(I)
    end
end
```
would store the on-grid value at each grid point of `itp` in the output `dest`.
Finally, courtesy of Julia's indexing rules, you can also use
```jl
fine = itp(linspace(1,10,1001), linspace(1,15,201))
```


## Control of interpolation algorithm

### BSplines

The interpolation type is described in terms of *degree*, *grid behavior* and, if necessary, *boundary conditions*. There are currently three degrees available: `Constant`, `Linear`, `Quadratic`,  and `Cubic` corresponding to B-splines of degree 0, 1, 2, and 3 respectively.

You also have to specify what *grid representation* you want. There are currently two choices: `OnGrid`, in which the supplied data points are assumed to lie *on* the boundaries of the interpolation interval, and `OnCell` in which the data points are assumed to lie on half-intervals between cell boundaries.

B-splines of quadratic or higher degree require solving an equation system to obtain the interpolation coefficients, and for that you must specify a *boundary condition* that is applied to close the system. The following boundary conditions are implemented: `Flat`, `Line` (alternatively, `Natural`), `Free`, `Periodic` and `Reflect`; their mathematical implications are described in detail in the pdf document under `/doc/latex`.

Some examples:
```jl
# Nearest-neighbor interpolation
itp = interpolate(a, BSpline(Constant()), OnCell())
v = itp(5.4)   # returns a[5]

# (Multi)linear interpolation
itp = interpolate(A, BSpline(Linear()), OnGrid())
v = itp(3.2, 4.1)  # returns 0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5])

# Quadratic interpolation with reflecting boundary conditions
# Quadratic is the lowest order that has continuous gradient
itp = interpolate(A, BSpline(Quadratic(Reflect())), OnCell())

# Linear interpolation in the first dimension, and no interpolation (just lookup) in the second
itp = interpolate(A, (BSpline(Linear()), NoInterp()), OnGrid())
v = itp(3.65, 5)  # returns  0.35*A[3,5] + 0.65*A[4,5]
```
There are more options available, for example:
```jl
# In-place interpolation
itp = interpolate!(A, BSpline(Quadratic(InPlace())), OnCell())
```
which destroys the input `A` but also does not need to allocate as much memory.

### Scaled BSplines

BSplines assume your data is uniformly spaced on the grid `1:N`, or its multidimensional equivalent. If you have data of the form `[f(x) for x in A]`, you need to tell Interpolations about the grid `A`. If `A` is not uniformly spaced, you must use gridded interpolation described below. However, if `A` is a collection of ranges or linspaces, you can use scaled BSplines. This is more efficient because the gridded algorithm does not exploit the uniform spacing. Scaled BSplines can also be used with any spline degree available for BSplines, while gridded interpolation does not currently support quadratic or cubic splines.

Some examples,
```jl
A_x = 1.:2.:40.
A = [log(x) for x in A_x]
itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())
sitp = scale(itp, A_x)
sitp(3.) # exactly log(3.)
sitp(3.5) # approximately log(3.5)
```

For multidimensional uniformly spaced grids
```jl
A_x1 = 1:.1:10
A_x2 = 1:.5:20
f(x1, x2) = log(x1+x2)
A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())
sitp = scale(itp, A_x1, A_x2)
sitp(5., 10.) # exactly log(5 + 10)
sitp(5.6, 7.1) # approximately log(5.6 + 7.1)
```
### Gridded interpolation

These use a very similar syntax to BSplines, with the major exception
being that one does not get to choose the grid representation (they
are all `OnGrid`). As such one must specify a set of coordinate arrays
defining the knots of the array.

In 1D
```jl
A = rand(20)
A_x = collect(1.0:2.0:40.0)
knots = (A_x,)
itp = interpolate(knots, A, Gridded(Linear()))
itp(2.0)
```

The spacing between adjacent samples need not be constant, you can use the syntax
```jl
itp = interpolate(knots, A, options...)
```
where `knots = (xknots, yknots, ...)` to specify the positions along
each axis at which the array `A` is sampled for arbitrary ("rectangular") samplings.

For example:
```jl
A = rand(8,20)
knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
itp = interpolate(knots, A, Gridded(Linear()))
itp(4,1.2)  # approximately A[2,6]
```
One may also mix modes, by specifying a mode vector in the form of an explicit tuple:
```jl
itp = interpolate(knots, A, (Gridded(Linear()),Gridded(Constant())))
```

Presently there are only three modes for gridded:
```jl
Gridded(Linear())
```
whereby a linear interpolation is applied between knots,
```jl
Gridded(Constant())
```
whereby nearest neighbor interpolation is used on the applied axis,
```jl
NoInterp
```
whereby the coordinate of the selected input vector MUST be located on a grid point. Requests for off grid
coordinates results in the throwing of an error.

## Parametric splines

Given a set a knots with coordinates `x(t)` and `y(t)`, a parametric spline `S(t) = (x(t),y(t))` parametrized by `t in [0,1]` can be constructed with the following code adapted from a [post](http://julia-programming-language.2336112.n4.nabble.com/Parametric-splines-td37794.html#a37818) by Tomas Lycken:

```julia
using Interpolations

t = 0:.1:1
x = sin.(2π*t)
y = cos.(2π*t)
A = hcat(x,y)

itp = scale(interpolate(A, (BSpline(Cubic(Natural())), NoInterp()), OnGrid()), t, 1:2)

tfine = 0:.01:1
xs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]
```

We can then plot the spline with:

```julia
using Plots

scatter(x, y, label="knots")
plot!(xs, ys, label="spline")
```
![parametric spline](doc/images/parametric_spline.png)

## Extrapolation

The call to `extrapolate` defines what happens if you try to index into the interpolation object with coordinates outside of `[1, size(data,d)]` in any dimension `d`. The implemented boundary conditions are `Throw`, `Flat`, `Linear`, `Periodic` and `Reflect`, with more options planned. `Periodic` and `Reflect` require that there is a method of `Base.mod` that can handle the indices used.

## Performance shootout

In the `perf` directory, you can find a script that tests
interpolation with several different packages.  We consider
interpolation in 1, 2, 3, and 4 dimensions, with orders 0
(`Constant`), 1 (`Linear`), and 2 (`Quadratic`).  Methods include
Interpolations `BSpline` (`IBSpline`) and `Gridded` (`IGridded`),
methods from the [Grid.jl](https://github.com/timholy/Grid.jl)
package, methods from the
[Dierckx.jl](https://github.com/kbarbary/Dierckx.jl) package, methods
from the
[GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
package (`GI`), methods from the
[ApproXD.jl](https://github.com/floswald/ApproXD.jl) package, and
methods from SciPy's `RegularGridInterpolator` accessed via `PyCall`
(`Py`).  All methods
are tested using an `Array` with approximately `10^6` elements, and
the interpolation task is simply to visit each grid point.

First, let's look at the two B-spline algorithms, `IBspline` and
`Grid`.  Here's a plot of the "construction time," the amount of time
it takes to initialize an interpolation object (smaller is better):

![construction](perf/constructionB.png)

The construction time is negligible until you get to second order
(quadratic); that's because quadratic is the lowest order requiring
the solution of tridiagonal systems upon construction.  The solvers
used by Interpolations are much faster than the approach taken in
Grid.

Now let's examine the interpolation performance.  Here we'll measure
"throughput", the number of interpolations performed per second
(larger is better):

![throughput](perf/rateB.png)

Once again, Interpolations wins on every test, by a factor that ranges
from 7 to 13.

Now let's look at the "gridded" methods that allow irregular spacing
along each axis.  For some of these, we compare interpolation performance in
both "vectorized" form `itp[xvector, yvector]` and in "scalar" form
`for y in yvector, x in xvector; val = itp[x,y]; end`.

First, construction time (smaller is better):

![construction](perf/constructionG.png)

Missing dots indicate cases that were not tested, or not supported by
the package.  (For construction, differences between "vec" and
"scalar" are just noise, since no interpolation is performed during
construction.)  The only package that takes appreciable construction
time is Dierckx.

And here's "throughput" (larger is better). To ensure we can see the
wide range of scales, here we use "square-root" scaling of the y-axis:

![throughput](perf/rateG.png)

For 1d, the "Dierckx scalar" and "GI" tests were interrupted because
they ran more than 20 seconds (far longer than any other test).  Both
performed much better in 2d, interestingly.  You can see that
Interpolations wins in every case, sometimes by a very large margin.

## Transitioning from Grid.jl

Instead of using
```julia
yi = InterpGrid(y, BCreflect, InterpQuadratic)
```
you should use
```julia
yi = interpolate(y, BSpline(Quadratic(Reflect())), OnCell())
```

In general, here are the closest mappings:

|  Grid             | Interpolations                             |
|:-----------------:|:------------------------------------------:|
| `InterpNearest`   | `Constant`                                 |
| `InterpLinear`    | `Linear`                                   |
| `InterpQuadratic` | `Quadratic`                                |
| `InterpCubic`     | `Cubic`                                    |
|                   |                                            |
| `BCnil`           | `extrapolate(itp, Interpolations.Throw())` |
| `BCnan`           | `extrapolate(itp, NaN)`                    |
| `BCna`            | `extrapolate(itp, NaN)`                    |
| `BCreflect`       | `interpolate` with `Reflect()`             |
| `BCperiodic`      | `interpolate` with `Periodic()`            |
| `BCnearest`       | `interpolate` with `Flat()`                |
| `BCfill`          | `extrapolate` with value                   |
|                   |                                            |
| odd orders        | `OnGrid()`                                 |
| even orders       | `OnCell()`                                 |


## Contributing

Work is very much in progress, but and help is always welcome. If you want to help out but don't know where to start, take a look at issue [#5 - our feature wishlist](https://github.com/JuliaMath/Interpolations.jl/issues/5) =) There is also some [developer documentation](doc/devdocs.md) that may help you understand how things work internally.

Contributions in any form are appreciated, but the best pull requests come with tests!
