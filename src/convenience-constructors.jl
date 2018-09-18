# convenience copnstructors for linear / cubic spline interpolations
# 1D version
LinearInterpolation(range::T, vs; extrapolation_bc = Throw()) where {T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Linear())), range), extrapolation_bc)
LinearInterpolation(range::T, vs; extrapolation_bc = Throw()) where {T <: AbstractArray} = extrapolate(interpolate((range, ), vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(range::T, vs; bc = Line(OnGrid()), extrapolation_bc = Throw()) where {T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), range), extrapolation_bc)

# multivariate versions
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Throw()) where {N,T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Linear())), ranges...), extrapolation_bc)
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Throw()) where {N,T <: AbstractArray} = extrapolate(interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(ranges::NTuple{N,T}, vs; bc = Line(OnGrid()), extrapolation_bc = Throw()) where {N,T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Cubic(bc))), ranges...), extrapolation_bc)
