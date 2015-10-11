export ScaledInterpolation

type ScaledInterpolation{T,N,ITPT,IT,GT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}
    itp::ITPT
    ranges::RT
end
ScaledInterpolation{T,ITPT,IT,GT,RT}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}, ranges::RT) =
    ScaledInterpolation{T,N,ITPT,IT,GT,RT}(itp, ranges)
"""
`scale(itp, xs, ys, ...)` scales an existing interpolation object to allow for indexing using other coordinate axes than unit ranges, by wrapping the interpolation object and transforming the indices from the provided axes onto unit ranges upon indexing.

The parameters `xs` etc must be either ranges or linspaces, and there must be one coordinate range/linspace for each dimension of the interpolation object.

For every `NoInterp` dimension of the interpolation object, the range must be exactly `1:size(itp, d)`.
"""
function scale{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ranges::Range...)
    length(ranges) == N || throw(ArgumentError("Must scale $N-dimensional interpolation object with exactly $N ranges (you used $(length(ranges)))"))
    for d in 1:N
        if iextract(IT,d) != NoInterp
            length(ranges[d]) == size(itp,d) || throw(ArgumentError("The length of the range in dimension $d ($(length(ranges[d]))) did not equal the size of the interpolation object in that direction ($(size(itp,d)))"))
        elseif ranges[d] != 1:size(itp,d)
            throw(ArgumentError("NoInterp dimension $d must be scaled with unit range 1:$(size(itp,d))"))
        end
    end

    ScaledInterpolation(T,N,itp,IT,GT,ranges)
end

@generated function getindex{T,N,ITPT,IT<:DimSpec}(sitp::ScaledInterpolation{T,N,ITPT,IT}, xs::Number...)
    length(xs) == N || throw(ArgumentError("Must index into $N-dimensional scaled interpolation object with exactly $N indices (you used $(length(xs)))"))
    interp_indices = map(i -> iextract(IT, i) != NoInterp ? :(coordlookup(sitp.ranges[$i], xs[$i])) : :(xs[$i]), 1:N)
    return :($(Expr(:meta,:inline)); getindex(sitp.itp, $(interp_indices...)))
end

getindex{T}(sitp::ScaledInterpolation{T,1}, x::Number, y::Int) = y == 1 ? sitp[x] : throw(BoundsError())

size(sitp::ScaledInterpolation, d) = size(sitp.itp, d)
lbound{T,N,ITPT,IT}(sitp::ScaledInterpolation{T,N,ITPT,IT,OnGrid}, d) = 1 <= d <= N ? sitp.ranges[d][1] : throw(BoundsError())
lbound{T,N,ITPT,IT}(sitp::ScaledInterpolation{T,N,ITPT,IT,OnCell}, d) = 1 <= d <= N ? sitp.ranges[d][1] - boundstep(sitp.ranges[d]) : throw(BoundsError())
ubound{T,N,ITPT,IT}(sitp::ScaledInterpolation{T,N,ITPT,IT,OnGrid}, d) = 1 <= d <= N ? sitp.ranges[d][end] : throw(BoundsError())
ubound{T,N,ITPT,IT}(sitp::ScaledInterpolation{T,N,ITPT,IT,OnCell}, d) = 1 <= d <= N ? sitp.ranges[d][end] + boundstep(sitp.ranges[d]) : throw(BoundsError())

boundstep(r::LinSpace) = ((r.stop - r.start) / r.divisor) / 2
boundstep(r::FloatRange) = r.step / 2
boundstep(r::StepRange) = r.step / 2
boundstep(r::UnitRange) = 1//2

"""
Returns *half* the width of one step of the range.

This function is used to calculate the upper and lower bounds of `OnCell` interpolation objects.
""" boundstep

coordlookup(r::LinSpace, x) = (r.divisor * x + r.stop - r.len * r.start) / (r.stop - r.start)
coordlookup(r::FloatRange, x) = (r.divisor * x - r.start) / r.step + one(eltype(r))
coordlookup(r::StepRange, x) = (x - r.start) / r.step + one(eltype(r))
coordlookup(r::UnitRange, x) = x - r.start + one(eltype(r))
coordlookup(i::Bool, r::Range, x) = i ? coordlookup(r, x) : convert(typeof(coordlookup(r,x)), x)

basetype{T,N,ITPT,IT,GT,RT}(::Type{ScaledInterpolation{T,N,ITPT,IT,GT,RT}}) = ITPT
basetype(sitp::ScaledInterpolation) = basetype(typeof(sitp))

gradient{T,N,ITPT,IT<:DimSpec}(sitp::ScaledInterpolation{T,N,ITPT,IT}, xs::Number...) = gradient!(Array(T,count_interp_dims(IT,N)), sitp, xs...)
@generated function gradient!{T,N,ITPT,IT}(g, sitp::ScaledInterpolation{T,N,ITPT,IT}, xs::Number...)
    ndims(g) == 1 || throw(DimensionMismatch("g must be a vector (but had $(ndims(g)) dimensions)"))
    length(xs) == N || throw(DimensionMismatch("Must index into $N-dimensional scaled interpolation object with exactly $N indices (you used $(length(xs)))"))

    interp_types = length(IT.parameters) == N ? IT.parameters : tuple([IT.parameters[1] for _ in 1:N]...)
    interp_dimens = map(it -> interp_types[it] != NoInterp, 1:N)
    interp_indices = map(i -> interp_dimens[i] ? :(coordlookup(sitp.ranges[$i], xs[$i])) : :(xs[$i]), 1:N)

    quote
        length(g) == $(count_interp_dims(IT, N)) || throw(ArgumentError(string("The length of the provided gradient vector (", length(g), ") did not match the number of interpolating dimensions (", $(count_interp_dims(IT, N)), ")")))
        gradient!(g, sitp.itp, $(interp_indices...))
        for i in eachindex(g)
            g[i] = rescale_gradient(sitp.ranges[i], g[i])
        end
        g
    end
end

rescale_gradient(r::LinSpace, g) = g * r.divisor / (r.stop - r.start)
rescale_gradient(r::FloatRange, g) = g * r.divisor / r.step
rescale_gradient(r::StepRange, g) = g / r.step
rescale_gradient(r::UnitRange, g) = g

"""
`rescale_gradient(r::Range)`

Implements the chain rule dy/dx = dy/du * du/dx for use when calculating gradients with scaled interpolation objects.
""" rescale_gradient
