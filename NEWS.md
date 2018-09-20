# v0.9.0

Breaking changes:

- `gradient` and `hessian` are no longer exported; use `Interpolations.gradient` and
  `Interpolations.hessian`.
- `interpolate` objects now check bounds, and throw an error if you try to evaluate them
  at locations beyond the edge of their interpolation domain; use `extrapolate` if you need out-of-bounds evaluation
- For quadratic and cubic interpolation, `interpolate!` now returns an object whose axes
  are narrowed by the amount of padding needed on the array edges. This preserves correspondence
  between input indices and output indices. See https://julialang.org/blog/2017/04/offset-arrays
  for more information.
- The parametrization of some types has changed; this does not affect users of the "exported"
  interface, but does break packages that performed manual construction of explicit types.

Changes with deprecation warnings:

- `itp[i...]` should be replaced with `itp(i...)`.
- `OnGrid` and `OnCell` should now be placed inside the boundary condition (e.g., `Flat(OnGrid())`),
  and should only be used for quadratic and cubic interpolation.
- the extrapolation boundary condition `Linear` was changed to `Line`, to be consistent
  with interpolation boundary conditions.

Advance notice of future changes:

- In future versions `itp[i...]` may be interpreted with reference to the parent array's
  indices rather than the knots supplied by the user (relevant for `scale` and `Gridded`).
  If you fix the existing deprecation warnings then you should be prepared for this change.
