
# Interpolations

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

```julia
using Pkg
Pkg.add("Interpolations")
```

## Example Usage
Create a grid `xs` and an array `A` of values to be interpolated
```julia
xs = 1:0.2:5
f(x) = log(x)
A = [f(x) for x in xs]
```
Create linear interpolation object without extrapolation
```julia
interp_linear = LinearInterpolation(xs, A)
interp_linear(3) # exactly log(3)
interp_linear(3.1) # approximately log(3.1)
interp_linear(0.9) # outside grid: error
```
Create linear interpolation object with extrapolation
```julia
interp_linear_extrap = LinearInterpolation(xs, A,extrapolation_bc=Line()) 
interp_linear_extrap(0.9) # outside grid: linear extrapolation
```
