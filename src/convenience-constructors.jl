# convenience copnstructors for constant / linear / cubic spline interpolations
# 1D version
ConstantInterpolation(range::AbstractRange, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Constant())), range), extrapolation_bc)
ConstantInterpolation(range::AbstractVector, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(interpolate((range, ), vs, Gridded(Constant())), extrapolation_bc)
LinearInterpolation(range::AbstractRange, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Linear())), range), extrapolation_bc)
LinearInterpolation(range::AbstractVector, vs::AbstractVector; extrapolation_bc = Throw()) =
    extrapolate(interpolate((range, ), vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(range::AbstractRange, vs::AbstractVector;
                         bc = Line(OnGrid()), extrapolation_bc = Throw()) =
    extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), range), extrapolation_bc)

# multivariate versions
ConstantInterpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Constant())), ranges...), extrapolation_bc)
ConstantInterpolation(ranges::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(interpolate(ranges, vs, Gridded(Constant())), extrapolation_bc)
LinearInterpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Linear())), ranges...), extrapolation_bc)
LinearInterpolation(ranges::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                    extrapolation_bc = Throw()) where {N,T} =
    extrapolate(interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(ranges::NTuple{N,AbstractRange}, vs::AbstractArray{T,N};
                         bc = Line(OnGrid()), extrapolation_bc = Throw()) where {N,T} =
    extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), ranges...), extrapolation_bc)

"""
    etp = LinearInterpolation(knots, A; extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, scheme), extrapolation_bc)`,
where `scheme` is either `BSpline(Linear())` or `Gridded(Linear())` depending on whether
`knots` are ranges or vectors.
"""
LinearInterpolation

"""
    etp = ConstantInterpolation(knots, A; extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, scheme), extrapolation_bc)`,
where `scheme` is either `BSpline(Constant())` or `Gridded(Constant())` depending on whether
`knots` are ranges or vectors.
"""
ConstantInterpolation

"""
    etp = CubicSplineInterpolation(knots, A; bc=Line(OnGrid()), extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, BSpline(Cubic(bc))), extrapolation_bc)`.
"""
CubicSplineInterpolation
