## Extrapolation

The call to `extrapolate` defines what happens if you try to index into the interpolation object with coordinates outside of its
bounds in any dimension. The implemented boundary conditions are `Throw`, `Flat`, `Linear`, `Periodic` and `Reflect`,
or you can pass a constant to be used as a "fill" value returned for any out-of-bounds evaluation.
`Periodic` and `Reflect` require that there is a method of `Base.mod` that can handle the indices used.

Examples:

```
itp = interpolate(1:7, BSpline(Linear()))
etpf = extrapolate(itp, Flat())   # gives 1 on the left edge and 7 on the right edge
etp0 = extrapolate(itp, 0)        # gives 0 everywhere outside [1,7]
```
