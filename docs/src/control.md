## Control of interpolation algorithm

### BSplines

The interpolation type is described in terms of *degree* and, if necessary, *boundary conditions*. There are currently four degrees available: `Constant`, `Linear`, `Quadratic`,  and `Cubic` corresponding to B-splines of degree 0, 1, 2, and 3 respectively.

B-splines of quadratic or higher degree require solving an equation system to obtain the interpolation coefficients, and for that you must specify a *boundary condition* that is applied to close the system. The following boundary conditions are implemented: `Flat`, `Line` (alternatively, `Natural`), `Free`, `Periodic` and `Reflect`; their mathematical implications are described in detail in the pdf document under `/doc/latex`.
When specifying these boundary conditions you also have to specify whether they apply at the edge grid point (`OnGrid()`)
or beyond the edge point halfway to the next (fictitious) grid point (`OnCell()`).

Some examples:
```julia
# Nearest-neighbor interpolation
itp = interpolate(a, BSpline(Constant()))
v = itp(5.4)   # returns a[5]

# (Multi)linear interpolation
itp = interpolate(A, BSpline(Linear()))
v = itp(3.2, 4.1)  # returns 0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5])

# Quadratic interpolation with reflecting boundary conditions
# Quadratic is the lowest order that has continuous gradient
itp = interpolate(A, BSpline(Quadratic(Reflect(OnCell()))))

# Linear interpolation in the first dimension, and no interpolation (just lookup) in the second
itp = interpolate(A, (BSpline(Linear()), NoInterp()))
v = itp(3.65, 5)  # returns  0.35*A[3,5] + 0.65*A[4,5]
```
There are more options available, for example:
```julia
# In-place interpolation
itp = interpolate!(A, BSpline(Quadratic(InPlace(OnCell()))))
```
which destroys the input `A` but also does not need to allocate as much memory.

### Scaled BSplines

BSplines assume your data is uniformly spaced on the grid `1:N`, or its multidimensional equivalent. If you have data of the form `[f(x) for x in A]`, you need to tell Interpolations about the grid `A`. If `A` is not uniformly spaced, you must use gridded interpolation described below. However, if `A` is a collection of ranges or linspaces, you can use scaled BSplines. This is more efficient because the gridded algorithm does not exploit the uniform spacing. Scaled BSplines can also be used with any spline degree available for BSplines, while gridded interpolation does not currently support quadratic or cubic splines.

Some examples,
```julia
A_x = 1.:2.:40.
A = [log(x) for x in A_x]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, A_x)
sitp(3.) # exactly log(3.)
sitp(3.5) # approximately log(3.5)
```

For multidimensional uniformly spaced grids
```julia
A_x1 = 1:.1:10
A_x2 = 1:.5:20
f(x1, x2) = log(x1+x2)
A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
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
```julia
A = rand(20)
A_x = collect(1.0:2.0:40.0)
knots = (A_x,)
itp = interpolate(knots, A, Gridded(Linear()))
itp(2.0)
```

The spacing between adjacent samples need not be constant, you can use the syntax
```julia
itp = interpolate(knots, A, options...)
```
where `knots = (xknots, yknots, ...)` to specify the positions along
each axis at which the array `A` is sampled for arbitrary ("rectangular") samplings.

For example:
```julia
A = rand(8,20)
knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
itp = interpolate(knots, A, Gridded(Linear()))
itp(4,1.2)  # approximately A[2,6]
```
One may also mix modes, by specifying a mode vector in the form of an explicit tuple:
```julia
itp = interpolate(knots, A, (Gridded(Linear()),Gridded(Constant())))
```

Presently there are only three modes for gridded:
```julia
Gridded(Linear())
```
whereby a linear interpolation is applied between knots,
```julia
Gridded(Constant())
```
whereby nearest neighbor interpolation is used on the applied axis,
```julia
NoInterp
```
whereby the coordinate of the selected input vector MUST be located on a grid point. Requests for off grid
coordinates results in the throwing of an error.

`missing` data will naturally propagate through the interpolation,
where some values will become missing. To avoid that, one can
filter out the missing data points and use a gridded interpolation.
For example:
```julia
x = 1:6
A = [i == 3 ? missing : i for i in x]
xf = [xi for (xi,a) in zip(x, A) if !ismissing(a)]
Af = [a for a in A if !ismissing(a)]
itp = interpolate((xf, ), Af, Gridded(Linear()))
```

## Parametric splines

Given a set a knots with coordinates `x(t)` and `y(t)`, a parametric spline `S(t) = (x(t),y(t))` parametrized by `t in [0,1]` can be constructed with the following code adapted from a [post](http://julia-programming-language.2336112.n4.nabble.com/Parametric-splines-td37794.html#a37818) by Tomas Lycken:

```julia
using Interpolations

t = 0:.1:1
x = sin.(2π*t)
y = cos.(2π*t)
A = hcat(x,y)

itp = scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)

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
