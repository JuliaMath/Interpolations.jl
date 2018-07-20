# convenience copnstructors for linear / cubic spline interpolations
# 1D version
LinearInterpolation(range::T, vs; extrapolation_bc = Interpolations.Throw()) where {T <: Range} = extrapolate(scale(interpolate(vs, BSpline(Linear()), OnGrid()), range), extrapolation_bc)
LinearInterpolation(range::T, vs; extrapolation_bc = Interpolations.Throw()) where {T <: AbstractArray} = extrapolate(interpolate((range, ), vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(range::T, vs; extrapolation_bc = Interpolations.Throw()) where {T <: Range} = extrapolate(scale(interpolate(vs, BSpline(Cubic(Line())), OnGrid()), range), extrapolation_bc)

# multivariate versions
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Interpolations.Throw()) where {N,T <: Range} = extrapolate(scale(interpolate(vs, BSpline(Linear()), OnGrid()), ranges...), extrapolation_bc)
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Interpolations.Throw()) where {N,T <: AbstractArray} = extrapolate(interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Interpolations.Throw()) where {N,T <: Range} = extrapolate(scale(interpolate(vs, BSpline(Cubic(Line())), OnGrid()), ranges...), extrapolation_bc)
