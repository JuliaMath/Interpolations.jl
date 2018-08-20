module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,
    scale,

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

using LinearAlgebra, SparseArrays
using StaticArrays, WoodburyMatrices, Ratios, AxisAlgorithms, OffsetArrays

using Base: @propagate_inbounds
import Base: convert, size, axes, promote_rule, ndims, eltype, checkbounds

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

"""
    BoundaryCondition

An abstract type with one of the following values (see the help for each for details):

- `Throw()`
- `Flat()`
- `Line()`
- `Free()`
- `Periodic()`
- `Reflect()`
- `InPlace()`
- `InPlaceQ()`
"""
abstract type BoundaryCondition <: Flag end
struct Throw <: BoundaryCondition end
struct Flat <: BoundaryCondition end
struct Line <: BoundaryCondition end
struct Free <: BoundaryCondition end
struct Periodic <: BoundaryCondition end
struct Reflect <: BoundaryCondition end
struct InPlace <: BoundaryCondition end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
struct InPlaceQ <: BoundaryCondition end
const Natural = Line

Base.IndexStyle(::Type{<:AbstractInterpolation}) = IndexCartesian()

size(itp::AbstractInterpolation) = map(length, itp.parentaxes)
size(exp::AbstractExtrapolation) = size(exp.itp)
axes(itp::AbstractInterpolation) = itp.parentaxes
axes(exp::AbstractExtrapolation) = axes(exp.itp)

bounds(itp::AbstractInterpolation) = map(twotuple, lbounds(itp), ubounds(itp))
twotuple(r::AbstractUnitRange) = (first(r), last(r))
twotuple(x, y) = (x, y)
bounds(itp::AbstractInterpolation, d) = bounds(itp)[d]
lbounds(itp::AbstractInterpolation) = _lbounds(itp.parentaxes, itpflag(itp), gridflag(itp))
ubounds(itp::AbstractInterpolation) = _ubounds(itp.parentaxes, itpflag(itp), gridflag(itp))
_lbounds(axs::NTuple{N,Any}, itp, gt) where N = ntuple(d->lbound(axs[d], iextract(itp,d), iextract(gt,d)), Val(N))
_ubounds(axs::NTuple{N,Any}, itp, gt) where N = ntuple(d->ubound(axs[d], iextract(itp,d), iextract(gt,d)), Val(N))

itptype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = IT
itptype(::Type{ITP}) where {ITP<:AbstractInterpolation} = itptype(supertype(ITP))
itptype(itp::AbstractInterpolation) = itptype(typeof(itp))
itpflag(itp::AbstractInterpolation) = itp.it
gridtype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = GT
gridtype(::Type{ITP}) where {ITP<:AbstractInterpolation} = gridtype(supertype(ITP))
gridtype(itp::AbstractInterpolation) = gridtype(typeof(itp))
gridflag(itp::AbstractInterpolation) = itp.gt
ndims(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = N
ndims(::Type{ITP}) where {ITP<:AbstractInterpolation} = ndims(supertype(ITP))
ndims(itp::AbstractInterpolation) = ndims(typeof(itp))
eltype(::Type{AbstractInterpolation{T,N,IT,GT}}) where {T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} = T
eltype(::Type{ITP}) where {ITP<:AbstractInterpolation} = eltype(supertype(ITP))
eltype(itp::AbstractInterpolation) = eltype(typeof(itp))

"""
    n = count_interp_dims(ITP)

Count the number of dimensions along which type `ITP` is interpolating.
`NoInterp` dimensions do not contribute to the sum.
"""
count_interp_dims(::Type{ITP}) where {ITP<:AbstractInterpolation} = count_interp_dims(itptype(ITP), ndims(ITP))
count_interp_dims(::Type{IT}, n) where {IT<:InterpolationType} = n * count_interp_dims(IT)
count_interp_dims(it::Type{IT}, n) where IT<:Tuple{Vararg{InterpolationType,N}} where N =
    _count_interp_dims(0, it...)
@inline _count_interp_dims(c, ::IT1, args...) where IT1 =
    _count_interp_dims(c + count_interp_dims(IT1), args...)
_count_interp_dims(c) = c

"""
    w = value_weights(degree, δx)

Compute the weights for interpolation of the value at an offset `δx` from the "base" position.
`degree` describes the interpolation scheme.

# Example

```jldoctest
julia> Interpolations.value_weights(Linear(), 0.2)
(0.8, 0.2)
```

This corresponds to the fact that linear interpolation at `x + 0.2` is `0.8*y[x] + 0.2*y[x+1]`.
"""
function value_weights end

"""
    w = gradient_weights(degree, δx)

Compute the weights for interpolation of the gradient at an offset `δx` from the "base" position.
`degree` describes the interpolation scheme.

# Example

```jldoctest
julia> Interpolations.gradient_weights(Linear(), 0.2)
(-1.0, 1.0)
```

This defines the gradient of a linear interpolation at 3.2 as `y[4] - y[3]`.
"""
function gradient_weights end

"""
    w = hessian_weights(degree, δx)

Compute the weights for interpolation of the hessian at an offset `δx` from the "base" position.
`degree` describes the interpolation scheme.

# Example

```jldoctest
julia> Interpolations.hessian_weights(Linear(), 0.2)
(0.0, 0.0)
```

Linear interpolation uses straight line segments, so the second derivative is zero.
"""
function hessian_weights end


gradient1(itp::AbstractInterpolation{T,1}, x) where {T} = gradient(itp, x)[1]
hessian1(itp::AbstractInterpolation{T,1}, x) where {T} = hessian(itp, x)[1]

### Supporting expansion of CartesianIndex

const UnexpandedIndexTypes = Union{Number, AbstractVector, CartesianIndex}
const ExpandedIndexTypes = Union{Number, AbstractVector}

Base.to_index(::AbstractInterpolation, x::Number) = x

# Commented out because you can't add methods to an abstract type.
# @inline function (itp::AbstractInterpolation)(x::Vararg{UnexpandedIndexTypes})
#     itp(to_indices(itp, x))
# end
function gradient(itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    gradient(itp, to_indices(itp, x)...)
end
function gradient!(dest, itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    gradient!(dest, itp, to_indices(itp, x)...)
end
function hessian(itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    hessian(itp, to_indices(itp, x)...)
end
function hessian!(dest, itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    hessian!(dest, itp, to_indices(itp, x)...)
end

# @inline function (itp::AbstractInterpolation)(x::Vararg{ExpandedIndexTypes})
#     itp.(Iterators.product(x...))
# end
function gradient(itp::AbstractInterpolation, x::Vararg{ExpandedIndexTypes})
    map(y->gradient(itp, y), Iterators.product(x...))
end
function hessian(itp::AbstractInterpolation, x::Vararg{ExpandedIndexTypes})
    map(y->hessian(itp, y), Iterators.product(x...))
end

# getindex is supported only for Integer indices (deliberately)
import Base: getindex
@propagate_inbounds getindex(itp::AbstractInterpolation{T,N}, i::Vararg{Integer,N}) where {T,N} = itp(i...)
@inline function getindex(itp::AbstractInterpolation{T,1}, i::Integer, j::Integer) where T
    @boundscheck (j == 1 || Base.throw_boundserror(itp, (i, j)))
    @inbounds itp(i)
end

# deprecate getindex for other numeric indices
@deprecate getindex(itp::AbstractInterpolation{T,N}, i::Vararg{Number,N}) where {T,N} itp(i...)

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
# include("gridded/gridded.jl")
# include("extrapolation/extrapolation.jl")
# include("scaling/scaling.jl")
include("utils.jl")
include("io.jl")
include("convenience-constructors.jl")

end # module
