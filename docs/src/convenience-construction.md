
## Convenience notation

For linear and cubic spline interpolations, `LinearInterpolation` and `CubicSplineInterpolation`
can be used to create interpolating and extrapolating objects handily:
```julia
f(x) = log(x)
xs = 1:0.2:5
A = [f(x) for x in xs]

# linear interpolation
interp_linear = LinearInterpolation(xs, A)
interp_linear(3) # exactly log(3)
interp_linear(3.1) # approximately log(3.1)

# cubic spline interpolation
interp_cubic = CubicSplineInterpolation(xs, A)
interp_cubic(3) # exactly log(3)
interp_cubic(3.1) # approximately log(3.1)
```
which support multidimensional data as well:
```julia
f(x,y) = log(x+y)
xs = 1:0.2:5
ys = 2:0.1:5
A = [f(x,y) for x in xs, y in ys]

# linear interpolation
interp_linear = LinearInterpolation((xs, ys), A)
interp_linear(3, 2) # exactly log(3 + 2)
interp_linear(3.1, 2.1) # approximately log(3.1 + 2.1)

# cubic spline interpolation
interp_cubic = CubicSplineInterpolation((xs, ys), A)
interp_cubic(3, 2) # exactly log(3 + 2)
interp_cubic(3.1, 2.1) # approximately log(3.1 + 2.1)
```
For extrapolation, i.e., when interpolation objects are evaluated in coordinates outside the range provided in constructors, the default option for a boundary condition is `Throw` so that they will return an error.
Interested users can specify boundary conditions by providing an extra parameter for `extrapolation_bc`:
```julia
f(x) = log(x)
xs = 1:0.2:5
A = [f(x) for x in xs]

# extrapolation with linear boundary conditions
extrap = LinearInterpolation(xs, A, extrapolation_bc = Line())

@test extrap(1 - 0.2) # ≈ f(1) - (f(1.2) - f(1))
@test extrap(5 + 0.2) # ≈ f(5) + (f(5) - f(4.8))
```
You can also use a "fill" value, which gets returned whenever you ask for out-of-range values:

```julia
extrap = LinearInterpolation(xs, A, extrapolation_bc = NaN)
@test isnan(extrap(5.2))
```

Irregular grids are supported as well; note that presently only `LinearInterpolation` supports irregular grids.
```julia
xs = [x^2 for x = 1:0.2:5]
A = [f(x) for x in xs]

# linear interpolation
interp_linear = LinearInterpolation(xs, A)
interp_linear(1) # exactly log(1)
interp_linear(1.05) # approximately log(1.05)
```
