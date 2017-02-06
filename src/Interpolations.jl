module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,
    scale,

    gradient!,
    gradient1,
    hessian!,
    hessian,
    hessian1,

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

using Compat
using WoodburyMatrices, Ratios, AxisAlgorithms

import Base: convert, size, getindex, gradient, promote_rule, ndims, eltype
if isdefined(:scaling) import Base.scaling end

abstract Flag
abstract InterpolationType <: Flag
immutable NoInterp <: InterpolationType end
abstract GridType <: Flag
immutable OnGrid <: GridType end
immutable OnCell <: GridType end

typealias DimSpec{T} Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} <: AbstractArray{T,N}
abstract AbstractInterpolationWrapper{T,N,ITPT,IT,GT} <: AbstractInterpolation{T,N,IT,GT}
abstract AbstractExtrapolation{T,N,ITPT,IT,GT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}

immutable Throw <: Flag end
immutable Flat <: Flag end
immutable Line <: Flag end
immutable Free <: Flag end
immutable Periodic <: Flag end
immutable Reflect <: Flag end
immutable InPlace <: Flag end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
immutable InPlaceQ <: Flag end
typealias Natural Line

# TODO: size might have to be faster?
size{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->size(itp,i), N)::NTuple{N,Int}
size(exp::AbstractExtrapolation, d) = size(exp.itp, d)
bounds{T,N}(itp::AbstractInterpolation{T,N}) = tuple(zip(lbounds(itp), ubounds(itp))...)
bounds{T,N}(itp::AbstractInterpolation{T,N}, d) = (lbound(itp,d),ubound(itp,d))
lbounds{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->lbound(itp,i), N)::NTuple{N,T}
ubounds{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->ubound(itp,i), N)::NTuple{N,T}
lbound{T,N}(itp::AbstractInterpolation{T,N}, d) = 1
ubound{T,N}(itp::AbstractInterpolation{T,N}, d) = size(itp, d)
itptype{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}}(::Type{AbstractInterpolation{T,N,IT,GT}}) = IT
itptype{ITP<:AbstractInterpolation}(::Type{ITP}) = itptype(supertype(ITP))
itptype(itp::AbstractInterpolation ) = itptype(typeof(itp))
gridtype{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}}(::Type{AbstractInterpolation{T,N,IT,GT}}) = GT
gridtype{ITP<:AbstractInterpolation}(::Type{ITP}) = gridtype(supertype(ITP))
gridtype(itp::AbstractInterpolation) = gridtype(typeof(itp))
ndims{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}}(::Type{AbstractInterpolation{T,N,IT,GT}}) = N
ndims{ITP<:AbstractInterpolation}(::Type{ITP}) = ndims(supertype(ITP))
ndims(itp::AbstractInterpolation) = ndims(typeof(itp))
eltype{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}}(::Type{AbstractInterpolation{T,N,IT,GT}}) = T
eltype{ITP<:AbstractInterpolation}(::Type{ITP}) = eltype(supertype(ITP))
eltype(itp::AbstractInterpolation) = eltype(typeof(itp))
count_interp_dims{T<:AbstractInterpolation}(::Type{T}, N) = N

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
include("gridded/gridded.jl")
include("extrapolation/extrapolation.jl")
include("scaling/scaling.jl")
include("utils.jl")

end # module