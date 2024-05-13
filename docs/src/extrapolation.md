## Extrapolation

The call to `extrapolate` defines what happens if you try to index into the
interpolation object with coordinates outside of its bounds in any dimension.
The implemented boundary conditions are `Throw`, `Flat`, `Line`, `Periodic` and
`Reflect`, or you can pass a constant to be used as a "fill" value returned for
any out-of-bounds evaluation.  `Periodic` and `Reflect` require that there is a
method of `Base.mod` that can handle the indices used.

Examples:

```julia
itp = interpolate(1:7, BSpline(Linear()))
etpf = extrapolate(itp, Flat())   # gives 1 on the left edge and 7 on the right edge
etp0 = extrapolate(itp, 0)        # gives 0 everywhere outside [1,7]
```

### Periodic extrapolation

For uniformly sampled periodic data, one can perform periodic extrapolation for all types of
B-Spline interpolations. By using the `Periodic(OnCell())` boundary condition in `interpolate`,
one does not need to include the periodic image of the starting sample point.

Examples:

```julia
f(x) = sin((x-3)*2pi/7 - 1)
A = Float64[f(x) for x in 1:7] # Does not include the periodic image

# Constant(Periodic())) is an alias for Constant{Nearest}(Periodic(OnCell()))
itp0 = interpolate(A, BSpline(Constant(Periodic())))
# Linear(Periodic())) is an alias for Linear(Periodic(OnCell()))
itp1 = interpolate(A, BSpline(Linear(Periodic())))
itp2 = interpolate(A, BSpline(Quadratic(Periodic(OnCell()))))
itp3 = interpolate(A, BSpline(Cubic(Periodic(OnCell()))))

etp0 = extrapolate(itp0, Periodic())
etp1 = extrapolate(itp1, Periodic())
etp2 = extrapolate(itp2, Periodic())
etp3 = extrapolate(itp3, Periodic())
```
