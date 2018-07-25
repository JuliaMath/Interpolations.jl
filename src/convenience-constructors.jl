# convenience copnstructors for linear / cubic spline interpolations
# 1D version
LinearInterpolation(range::T, vs; extrapolation_bc = Interpolations.Throw()) where {T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Linear()), OnGrid()), range), extrapolation_bc)
LinearInterpolation(range::T, vs; extrapolation_bc = Interpolations.Throw()) where {T <: AbstractArray} = extrapolate(interpolate((range, ), vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(range::T, vs; bc = Interpolations.Line(), extrapolation_bc = Interpolations.Throw()) where {T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Cubic(bc)), OnGrid()), range), extrapolation_bc)

# multivariate versions
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Interpolations.Throw()) where {N,T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Linear()), OnGrid()), ranges...), extrapolation_bc)
LinearInterpolation(ranges::NTuple{N,T}, vs; extrapolation_bc = Interpolations.Throw()) where {N,T <: AbstractArray} = extrapolate(interpolate(ranges, vs, Gridded(Linear())), extrapolation_bc)
CubicSplineInterpolation(ranges::NTuple{N,T}, vs; bc = Interpolations.Line(), extrapolation_bc = Interpolations.Throw()) where {N,T <: AbstractRange} = extrapolate(scale(interpolate(vs, BSpline(Cubic(bc)), OnGrid()), ranges...), extrapolation_bc)
