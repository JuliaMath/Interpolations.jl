# convenience copnstructors for constant / linear / cubic spline interpolations
# 1D version
constant_interpolation(range::AbstractRange, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Constant())), range), extrapolation_bc)
constant_interpolation(range::AbstractVector, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(interpolate((range, ), vs, Gridded(Constant())), extrapolation_bc)
linear_interpolation(range::AbstractRange, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Linear())), range), extrapolation_bc)
linear_interpolation(range::AbstractVector, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(interpolate((range, ), vs, Gridded(Linear())), extrapolation_bc)
cubic_spline_interpolation(range::AbstractRange, vs::AbstractVector;
                         bc = Line(OnGrid()), extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), range), extrapolation_bc)

# multivariate versions
constant_interpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Constant())), ranges...), extrapolation_bc)
constant_interpolation(ranges::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(interpolate(ranges, vs, Gridded(Constant())), extrapolation_bc)
linear_interpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Linear())), ranges...), extrapolation_bc)
linear_interpolation(ranges::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
cubic_spline_interpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                         bc = Line(OnGrid()), extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), ranges...), extrapolation_bc)

"""
    etp = linear_interpolation(knots, A; extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, scheme), extrapolation_bc)`,
where `scheme` is either `BSpline(Linear())` or `Gridded(Linear())` depending on whether
`knots` are ranges or vectors.
"""
linear_interpolation

"""
    etp = constant_interpolation(knots, A; extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, scheme), extrapolation_bc)`,
where `scheme` is either `BSpline(Constant())` or `Gridded(Constant())` depending on whether
`knots` are ranges or vectors.
"""
constant_interpolation

"""
    etp = cubic_spline_interpolation(knots, A; bc=Line(OnGrid()), extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, BSpline(Cubic(bc))), extrapolation_bc)`.
"""
cubic_spline_interpolation
