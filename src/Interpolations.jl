module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,
    scale,

    gradient!,

    AbstractInterpolation,
    AbstractExtrapolation,

    OnCell,
    OnGrid,

    Flat,
    Line,
    Free,
    Periodic,
    Reflect,
    Natural,
    InPlace,
    InPlaceQ

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl
    # scaling/scaling.jl

using WoodburyMatrices, Ratios, AxisAlgorithms

import Base: convert, size, getindex, gradient, scale, promote_rule

abstract InterpolationType
immutable NoInterp <: InterpolationType end
abstract GridType
immutable OnGrid <: GridType end
immutable OnCell <: GridType end

typealias DimSpec{T} Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} <: AbstractArray{T,N}
abstract AbstractInterpolationWrapper{T,N,ITPT,IT,GT} <: AbstractInterpolation{T,N,IT,GT}
abstract AbstractExtrapolation{T,N,ITPT,IT,GT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}

abstract BoundaryCondition
immutable Flat <: BoundaryCondition end
immutable Line <: BoundaryCondition end
immutable Free <: BoundaryCondition end
immutable Periodic <: BoundaryCondition end
immutable Reflect <: BoundaryCondition end
immutable InPlace <: BoundaryCondition end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
immutable InPlaceQ <: BoundaryCondition end
typealias Natural Line

# TODO: size might have to be faster?
size{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->size(itp,i), N)::NTuple{N,Int}
size(exp::AbstractExtrapolation, d) = size(exp.itp, d)
bounds{T,N}(itp::AbstractInterpolation{T,N}) = tuple(zip(lbounds(itp), ubounds(itp))...)
bounds{T,N}(itp::AbstractInterpolation{T,N}, d) = (lbound(itp,d),ubound(itp,d))
lbounds{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->lbound(itp,i), N)::NTuple{N,T}
ubounds{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->ubound(itp,i), N)::NTuple{N,T}
lbound{T,N}(itp::AbstractInterpolation{T,N}, d) = convert(T, 1)
ubound{T,N}(itp::AbstractInterpolation{T,N}, d) = convert(T, size(itp, d))
itptype{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}) = IT
gridtype{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}) = GT
count_interp_dims{T<:AbstractInterpolation}(::Type{T}, N) = N

@generated function gradient{T,N}(itp::AbstractInterpolation{T,N}, xs...)
    n = count_interp_dims(itp, N)
    Tg = promote_type(T, [x <: AbstractArray ? eltype(x) : x for x in xs]...)
    xargs = [:(xs[$d]) for d in 1:length(xs)]
    :(gradient!(Array($Tg,$n), itp, $(xargs...)))
end

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
include("gridded/gridded.jl")
include("extrapolation/extrapolation.jl")
include("scaling/scaling.jl")

end # module
