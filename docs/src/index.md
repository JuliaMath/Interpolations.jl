
# Interpolations

[![version](https://juliahub.com/docs/Interpolations/version.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx)
[![pkgeval](https://juliahub.com/docs/Interpolations/pkgeval.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx)
[![Build Status](https://github.com/JuliaMath/Interpolations.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/Interpolations.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![deps](https://juliahub.com/docs/Interpolations/deps.svg)](https://juliahub.com/ui/Packages/Interpolations/VpKVx?t=2)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://juliamath.github.io/Interpolations.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliamath.github.io/Interpolations.jl/latest)

The package [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl/)
implements a variety of interpolation schemes for the
Julia language.  It has the goals of ease-of-use, broad algorithmic
support, and exceptional performance.

Currently this package supports
[B-splines](https://en.wikipedia.org/wiki/B-spline) and
irregular grids.  The API has been designed with
intent to support more options. Pull-requests are more than welcome!
It should be noted that the API may continue to evolve over time.

There are many other interpolation packages implemented in Julia.
For a listing, see [Other Interpolation Packages](@ref).

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

## Installation

Interpolations.jl can be installed via the following invocation
since it is a registered Julia package.

```julia-repl
julia> using Pkg
julia> Pkg.add("Interpolations")
```

## Example Usage
Create a grid `xs` and an array `A` of values to be interpolated
```@example lerp
using Interpolations
xs = 1:0.2:5
A = log.(xs)
nothing # hide
```
Create linear interpolation object without extrapolation
```@repl lerp
interp_linear = linear_interpolation(xs, A);
interp_linear(3) # exactly log(3)
interp_linear(3.1) # approximately log(3.1)
interp_linear(0.9) # outside grid: error
```
Create linear interpolation object with extrapolation
```@repl lerp
interp_linear_extrap = linear_interpolation(xs, A, extrapolation_bc=Line());
interp_linear_extrap(0.9) # outside grid: linear extrapolation
```


## Performant Example Usage

The above use of `linear_interpolation` is actually a short hand for a
composition of `interpolate`, `scale`, and `extrapolate`. You may not need all
of the the scaling and extrapolation features.

```julia
interp_linear = extrapolate(scale(interpolate(A, BSpline(Linear())), xs), Line())
```

If we know we do not need the extrapolation portion, we can use the following.

```julia
scaled_itp = scale(interpolate(A, BSpline(Linear())), xs)
```

We can also remove the scaling for further performance if integer valued knots
and regular grids are sufficient.

```julia
itp = interpolate(A, BSpline(Linear()))
```

Removing the scaling or extrapolation will help accelerate interpolation by
removing unneeded operations and branches. This can permit the use of advanced
processor Single Instruction/Multiple Data (SIMD) capabilities.

## Regular Grids

Interpolations.jl is optimized for use with regular grids with uniform spacing.
The highest performance is achieved when the knots are an `AbstractUnitRange`
such as `2:5` or `Base.OneTo(9)`. The default case if no knots are specified is to
assign the knots as a `UnitRange` starting at `1`.

## Scaling

If the knots are not unit spaced or start at a distinct value other than `1`,
then the `scale` function can be used. While this increases the flexibility of
the interpolation, some performance penalty is acquired.
See [Scaled BSplines](@ref) for further information.

## Irregular Grids

If the knots are irregularly spaced, then the ranges between knots will have to
be scaled as in the `Gridded` interpolation type. See [Gridded interpolation](@ref)
for additional details.

## Points outside the knots

For points not between knots, extrapolation can be used. This introduces a
branch into the code that checks whether the point to be queried is inside or
outside of the knots. This branch can inhibit the use of vectorized SIMD
computation, resulting in a reduction of performance. See [Extrapolation](@ref).
