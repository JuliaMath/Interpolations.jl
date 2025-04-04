# Interpolations.jl

[![version](https://juliahub.com/docs/Interpolations/version.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx)
[![pkgeval](https://juliahub.com/docs/Interpolations/pkgeval.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx)
[![Build Status](https://github.com/JuliaMath/Interpolations.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/Interpolations.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![deps](https://juliahub.com/docs/Interpolations/deps.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx?t=2)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://juliamath.github.io/Interpolations.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliamath.github.io/Interpolations.jl/latest)

This package implements a variety of interpolation schemes for the
Julia language.  It has the goals of ease-of-use, broad algorithmic
support, and exceptional performance.

Currently this package supports
[B-splines](https://en.wikipedia.org/wiki/B-spline) and also
irregular grids.  The API has been designed with
intent to support more options. Initial support for Lanczos
interpolation was recently added. Pull-requests are more than welcome!
It should be noted that the API may continue to evolve over time.

At the bottom of this page, you can find a "performance shootout"
among these methods (as well as SciPy's `RegularGridInterpolator`).

## Installation

Interpolations.jl can be installed via the following invocation
since it is a registered Julia package.

```julia
using Pkg
Pkg.add("Interpolations")
```

## Example Usage
Create a grid `xs` and an array `A` of values to be interpolated
```julia
xs = 1:0.2:5
A = log.(xs)
```
Create linear interpolation object without extrapolation
```julia
interp_linear = linear_interpolation(xs, A)
interp_linear(3) # exactly log(3)
interp_linear(3.1) # approximately log(3.1)
interp_linear(0.9) # outside grid: error
```
Create linear interpolation object with extrapolation
```julia
interp_linear_extrap = linear_interpolation(xs, A,extrapolation_bc=Line()) 
interp_linear_extrap(0.9) # outside grid: linear extrapolation
```

### Other Examples

More examples, such as plotting and cubic interpolation, can be found at the [convenience constructions](docs/src/convenience-construction.md#example-with-plotsjl) documentation.

![interpolation plot example](docs/src/assets/plotsjl_interpolation_example.png)

## Other Interpolation Packages

Other interpolation packages for Julia include:

- [ApproXD.jl](https://github.com/floswald/ApproXD.jl) implements B-spline and linear interpolation in Julia.
- [BarycentricInterpolation.jl](https://github.com/dawbarton/BarycentricInterpolation.jl) implements the Barycentric formula for polynomial interpolation on equispaced points and Chebyshev points of the first and second kind.
- [BasicInterpolators.jl](https://github.com/markmbaum/BasicInterpolators.jl) provides a collection of common interpolation recipes for basic applications.
- [BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) offers tools for B-spline based Galerkin and collocation methods, including for interpolation and approximation.
- [Curves.jl](https://github.com/lungben/Curves.jl) supports log-interpolation via immutable `Curve` objects.
- [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) is a library for performing interpolations of one-dimensional data.
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl) is a wrapper for the dierckx Fortran library, which also underlies `scipy.interpolate`.
- [DIVAnd.jl](https://github.com/gher-ulg/DIVAnd.jl) for N-dimensional smoothing interpolation. 
- [FastChebInterp.jl](https://github.com/stevengj/FastChebInterp.jl) does fast multidimensional Chebyshev interpolation on a hypercube using separable grid of interpolation points.
- [FEMBasis.jl](https://github.com/JuliaFEM/FEMBasis.jl) contains interpolation routines for standard finite element function spaces.
- [FineShift.jl](https://github.com/emmt/FineShift.jl) does fast sub-sample shifting of multidimensional arrays.
- [FourierTools.jl](https://github.com/bionanoimaging/FourierTools.jl) includes sinc interpolation for up and down sampling.
- [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) provides interpolation and simulation methods over complex 2D and 3D meshes.
- [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl) performs multivariate interpolation on a rectilinear grid.
- [InterpolationKernels.jl](https://github.com/emmt/InterpolationKernels.jl) provides a library of interpolation kernels.
- [KernelInterpolation.jl](https://github.com/JoshuaLampert/KernelInterpolation.jl) implements scattered data interpolations in arbitrary dimensions by radial basis functions with support for solving linear partial differential equations.
- [KissSmoothing.jl](https://github.com/francescoalemanno/KissSmoothing.jl) implements denoising and a Radial Basis Function estimation procedure.
- [LinearInterpolations.jl](https://github.com/jw3126/LinearInterpolations.jl) allows for interpolation using weighted averages allowing probability distributions, rotations, and other Lie groups to be interpolated.
- [LinearInterpolators.jl](https://github.com/emmt/LinearInterpolators.jl) provides linear interpolation methods for Julia based on InterpolationKernels.jl, above.
- [LocalFunctionApproximation.jl](https://github.com/sisl/LocalFunctionApproximation.jl) provides local function approximators that interpolates a scalar-valued function across a vector space.
- [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) provides natural neighbour interpolation methods for scattered two-dimensional point sets, with support for derivative generation.
- [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) for monotonic interpolation.
- [PiecewiseLinearApprox.jl](https://github.com/RJDennis/PiecewiseLinearApprox.jl) performs piecewise linear interpolation over an arbitrary number of dimensions.
- [ScatteredInterpolation.jl](https://github.com/eljungsk/ScatteredInterpolation.jl) interpolates scattered data in arbitrary dimensions.

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

If you would like to list a registered package that is related to interpolation, please create a Github issue.


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

## Development Status

This package is being maintained but not actively developed. Maintenance is
focused on fixing bugs and issues with the current code base. New features are
welcome via pull requests and will be reviewed and released in a timely fashion.

If you would like to become involved in maintenance or active development of
the package please feel free to get in touch via a Github issue.

This package follows semantic version in that documented features should not
break without changing the minor version.

See the [news](NEWS.md) for details on how to update between breaking releases,
indicated by changes in minor versions.

## Contributing

Work is very much in progress, but help is always welcome. There is some [developer documentation](http://juliamath.github.io/Interpolations.jl/latest/devdocs/) that may help you understand how things work internally.

Contributions in any form are appreciated, but the best pull requests come with tests!
