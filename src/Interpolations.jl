VERSION < v"0.7.0-beta2.199" && __precompile__()

module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,
    scale,

    gradient,
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
    InPlaceQ,

    LinearInterpolation,
    CubicSplineInterpolation

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl
    # scaling/scaling.jl

using Compat
using Compat.LinearAlgebra, Compat.SparseArrays
using WoodburyMatrices, Ratios, AxisAlgorithms, OffsetArrays

import Base: convert, size, getindex, promote_rule,
             ndims, eltype, checkbounds

@static if VERSION < v"0.7.0-DEV.3449"
    import Base: gradient
end

import Compat: axes

abstract type Flag end
abstract type InterpolationType <: Flag end
struct NoInterp <: InterpolationType end
abstract type GridType <: Flag end
struct OnGrid <: GridType end
struct OnCell <: GridType end

const DimSpec{T} = Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract type AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} <: AbstractArray{T,N} end
abstract type AbstractInterpolationWrapper{T,N,ITPT,IT,GT} <: AbstractInterpolation{T,N,IT,GT} end
abstract type AbstractExtrapolation{T,N,ITPT,IT,GT} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT} end

struct Throw <: Flag end
struct Flat <: Flag end
struct Line <: Flag end
struct Free <: Flag end
struct Periodic <: Flag end
struct Reflect <: Flag end
struct InPlace <: Flag end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
struct InPlaceQ <: Flag end
const Natural = Line

@generated size(itp::AbstractInterpolation{T,N}) where {T, N} = Expr(:tuple, [:(size(itp, $i)) for i in 1:N]...)
size(exp::AbstractExtrapolation, d) = size(exp.itp, d)
bounds(itp::AbstractInterpolation{T,N}) where {T,N} = tuple(zip(lbounds(itp), ubounds(itp))...)
bounds(itp::AbstractInterpolation{T,N}, d) where {T,N} = (lbound(itp,d),ubound(itp,d))
@generated lbounds(itp::AbstractInterpolation{T,N}) where {T,N} = Expr(:tuple, [:(lbound(itp, $i)) for i in 1:N]...)
@generated ubounds(itp::AbstractInterpolation{T,N}) where {T,N} = Expr(:tuple, [:(ubound(itp, $i)) for i in 1:N]...)
lbound(itp::AbstractInterpolation{T,N}, d) where {T,N} = 1
ubound(itp::AbstractInterpolation{T,N}, d) where {T,N} = size(itp, d)
itptype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = IT
itptype(::Type{ITP}) where {ITP<:AbstractInterpolation} = itptype(supertype(ITP))
itptype(itp::AbstractInterpolation ) = itptype(typeof(itp))
gridtype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = GT
gridtype(::Type{ITP}) where {ITP<:AbstractInterpolation} = gridtype(supertype(ITP))
gridtype(itp::AbstractInterpolation) = gridtype(typeof(itp))
ndims(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = N
ndims(::Type{ITP}) where {ITP<:AbstractInterpolation} = ndims(supertype(ITP))
ndims(itp::AbstractInterpolation) = ndims(typeof(itp))
eltype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = T
eltype(::Type{ITP}) where {ITP<:AbstractInterpolation} = eltype(supertype(ITP))
eltype(itp::AbstractInterpolation) = eltype(typeof(itp))
count_interp_dims(::Type{T}, N) where {T<:AbstractInterpolation} = N

# Generic indexing methods (fixes #183)
Base.to_index(::AbstractInterpolation, i::Real) = i

@inline function Base._getindex(::IndexCartesian, A::AbstractInterpolation{T,N}, I::Vararg{Int,N}) where {T,N} # ambiguity resolution
    @inbounds r = getindex(A, I...)
    r
end
@inline function Base._getindex(::IndexCartesian, A::AbstractInterpolation{T,N}, I::Vararg{Real,N}) where {T,N}
    @inbounds r = getindex(A, I...)
    r
end

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
include("gridded/gridded.jl")
include("extrapolation/extrapolation.jl")
include("scaling/scaling.jl")
include("utils.jl")
include("io.jl")
include("convenience-constructors.jl")

end # module
