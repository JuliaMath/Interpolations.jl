# convenience constructors for constant / linear / cubic spline interpolations
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

# convolution interpolation
# convolution_interpolation(knots::Union{AbstractVector,NTuple{N,AbstractVector}}, values::AbstractArray{T}; 
#     extrapolation_bc = Throw()) where {T,N} = extrapolate(scale(ConvolutionInterpolation(knots isa AbstractVector ? (knots,) : knots, values), knots...), extrapolation_bc)
# convolution_interpolation(ranges::NTuple{N,AbstractVector}, vs::AbstractArray{T}; 
#     extrapolation_bc = Throw()) where {T,N} = ConvolutionInterpolation(ranges isa AbstractVector ? (ranges,) : ranges, vs) # extrapolate(, extrapolation_bc)
function convolution_interpolation(knots::Union{AbstractVector,NTuple{N,AbstractVector}}, values::AbstractArray{T}; 
    extrapolation_bc = Throw()) where {T,N}
    if knots isa AbstractVector
        knots = (knots,)
    end
    itp = ConvolutionInterpolation(knots, values)
    return extrapolate(itp, extrapolation_bc)
end
"""
    etp = linear_interpolation(knots, A; extrapolation_bc=Throw())

A shorthand for one of the following.
* `extrapolate(scale(interpolate(A, BSpline(Linear())), knots), extrapolation_bc)`
* `extrapolate(interpolate(knots, A, Gridded(Linear())), extrapolation_bc)`,
depending on whether `knots` are ranges or vectors.

Consider using `interpolate`, `scale`, or `extrapolate` explicitly as needed
rather than using this convenience constructor. Performance will improve
without scaling or extrapolation.
"""
linear_interpolation

"""
    etp = constant_interpolation(knots, A; extrapolation_bc=Throw())

A shorthand for `extrapolate(interpolate(knots, A, scheme), extrapolation_bc)`,
where `scheme` is either `BSpline(Constant())` or `Gridded(Constant())` depending on whether
`knots` are ranges or vectors.

Consider using `interpolate` or `extrapolate` explicitly as needed
rather than using this convenience constructor. Performance will improve
without extrapolation.
"""
constant_interpolation

"""
    etp = cubic_spline_interpolation(knots, A; bc=Line(OnGrid()), extrapolation_bc=Throw())

A shorthand for `extrapolate(scale(interpolate(A, BSpline(Cubic(bc))), knots), extrapolation_bc)`.

Consider using `interpolate`, `scale`, or `extrapolate` explicitly as needed
rather than using this convenience constructor. Performance will improve
without scaling or extrapolation.
"""
cubic_spline_interpolation
"""
    etp = convolution_interpolation(knots, values; extrapolation_bc=Throw())

A shorthand for `extrapolate(ConvolutionInterpolation(knots, values), extrapolation_bc)`.
"""
convolution_interpolation