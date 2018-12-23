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
    Throw,

    LinearInterpolation,
    CubicSplineInterpolation

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl
    # monotonic/monotonic.jl
    # scaling/scaling.jl

using LinearAlgebra, SparseArrays
using StaticArrays, WoodburyMatrices, Ratios, AxisAlgorithms, OffsetArrays

using Base: @propagate_inbounds
import Base: convert, size, axes, promote_rule, ndims, eltype, checkbounds, axes1

abstract type Flag end
abstract type InterpolationType <: Flag end
struct NoInterp <: InterpolationType end
abstract type GridType <: Flag end
struct OnGrid <: GridType end
struct OnCell <: GridType end

const DimSpec{T} = Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract type AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType}} <: AbstractArray{T,N} end
abstract type AbstractInterpolationWrapper{T,N,ITPT,IT} <: AbstractInterpolation{T,N,IT} end
abstract type AbstractExtrapolation{T,N,ITPT,IT} <: AbstractInterpolationWrapper{T,N,ITPT,IT} end

"""
    BoundaryCondition

An abstract type with one of the following values (see the help for each for details):

- `Throw(gt)`
- `Flat(gt)`
- `Line(gt)`
- `Free(gt)`
- `Periodic(gt)`
- `Reflect(gt)`
- `InPlace(gt)`
- `InPlaceQ(gt)`

where `gt` is the grid type, e.g., `OnGrid()` or `OnCell()`. `OnGrid` means that the boundary
condition "activates" at the first and/or last integer location within the interpolation region,
`OnCell` means the interpolation extends a half-integer beyond the edge before
activating the boundary condition.
"""
abstract type BoundaryCondition <: Flag end
# Put the gridtype into the boundary condition, since that's all it affects (see issue #228)
# Nothing is used for extrapolation
struct Throw{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct Flat{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct Line{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct Free{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct Periodic{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct Reflect{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
struct InPlace{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
struct InPlaceQ{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
const Natural = Line

(::Type{BC})() where BC<:BoundaryCondition = BC(nothing)
function Base.show(io::IO, bc::BoundaryCondition)
    print(io, nameof(typeof(bc)), '(')
    bc.gt === nothing || show(io, bc.gt)
    print(io, ')')
end


Base.IndexStyle(::Type{<:AbstractInterpolation}) = IndexCartesian()

size(exp::AbstractExtrapolation) = size(exp.itp)
axes(exp::AbstractExtrapolation) = axes(exp.itp)

twotuple(r::AbstractUnitRange) = (first(r), last(r))
twotuple(x, y) = (x, y)
bounds(itp::AbstractInterpolation) = map(twotuple, lbounds(itp), ubounds(itp))
bounds(itp::AbstractInterpolation, d) = bounds(itp)[d]

itptype(::Type{AbstractInterpolation{T,N,IT}}) where {T,N,IT<:DimSpec{InterpolationType}} = IT
itptype(::Type{ITP}) where {ITP<:AbstractInterpolation} = itptype(supertype(ITP))
itptype(itp::AbstractInterpolation) = itptype(typeof(itp))
gridtype(::Type{AbstractInterpolation{T,N,IT}}) where {T,N,IT<:DimSpec{InterpolationType}} = GT
gridtype(::Type{ITP}) where {ITP<:AbstractInterpolation} = gridtype(supertype(ITP))
gridtype(itp::AbstractInterpolation) = gridtype(typeof(itp))
ndims(::Type{AbstractInterpolation{T,N,IT}}) where {T,N,IT<:DimSpec{InterpolationType}} = N
ndims(::Type{ITP}) where {ITP<:AbstractInterpolation} = ndims(supertype(ITP))
ndims(itp::AbstractInterpolation) = ndims(typeof(itp))
eltype(::Type{AbstractInterpolation{T,N,IT}}) where {T,N,IT<:DimSpec{InterpolationType}} = T
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
    BoundsCheckStyle(itp)

A trait to determine dispatch of bounds-checking for `itp`.
Can return `NeedsCheck()`, in which case bounds-checking is performed, or `CheckWillPass()`
in which case the check will return `true`.
"""
abstract type BoundsCheckStyle end
struct NeedsCheck <: BoundsCheckStyle end
struct CheckWillPass <: BoundsCheckStyle end

BoundsCheckStyle(itp) = NeedsCheck()

"""
    wi = WeightedIndex(indexes, weights)

Construct a weighted index `wi`, which can be thought of as a generalization of an
ordinary array index to the context of interpolation.
For an ordinary vector `a`, `a[i]` extracts the element at index `i`.
When interpolating, one is typically interested in a range of indexes and the output is
some weighted combination of array values at these indexes.
For example, for linear interpolation at a location `x` between integers `i` and `i+1`, we have

    ret = (1-f)*a[i] + f*a[i+1]

where `f = x-i` lies between 0 and 1. This can be represented as `a[wi]`, where

    wi = WeightedIndex(i:i+1, (1-f, f))

i.e.,

    ret = sum(a[indexes] .* weights)

Linear interpolation thus constructs weighted indices using a 2-tuple for `weights` and
a length-2 `indexes` range.
Higher-order interpolation would involve more positions and weights (e.g., 3-tuples for
quadratic interpolation, 4-tuples for cubic).

In multiple dimensions, separable interpolation schemes are implemented in terms
of multiple weighted indices, accessing `A[wi1, wi2, ...]` where each `wi` is the
`WeightedIndex` along the corresponding dimension.

For value interpolation, `weights` will typically sum to 1.
However, for gradient and Hessian computation this will not necessarily be true.
For example, the gradient of one-dimensional linear interpolation can be represented as

    gwi = WeightedIndex(i:i+1, (-1, 1))
    g1 = a[gwi]

For a three-dimensional array `A`, one might compute `∂A/∂x₂` (the second component
of the gradient) as `A[wi1, gwi2, wi3]`, where `wi1` and `wi3` are "value" weights
and `gwi2` "gradient" weights.

`indexes` may be supplied as a range or as a tuple of the same length as `weights`.
The latter is applicable, e.g., for periodic boundary conditions.
"""
abstract type WeightedIndex{L,W} end

# Type to use when array locations are adjacent. This may offer more opportunities
# for compiler optimizations (e.g., SIMD).
struct WeightedAdjIndex{L,W} <: WeightedIndex{L,W}
    istart::Int
    weights::NTuple{L,W}
end
# Type to use with non-adjacent locations. E.g., periodic boundary conditions.
struct WeightedArbIndex{L,W} <: WeightedIndex{L,W}
    indexes::NTuple{L,Int}
    weights::NTuple{L,W}
end

function WeightedIndex(indexes::AbstractUnitRange{<:Integer}, weights::NTuple{L,Any}) where L
    @noinline mismatch(indexes, weights) = throw(ArgumentError("the length of indexes must match weights, got $indexes vs $weights"))
    length(indexes) == L || mismatch(indexes, weights)
    WeightedAdjIndex(first(indexes), promote(weights...))
end
WeightedIndex(istart::Integer, weights::NTuple{L,Any}) where L =
    WeightedAdjIndex(istart, promote(weights...))
WeightedIndex(indexes::NTuple{L,Integer}, weights::NTuple{L,Any}) where L =
    WeightedArbIndex(indexes, promote(weights...))

weights(wi::WeightedIndex) = wi.weights
indexes(wi::WeightedAdjIndex) = wi.istart
indexes(wi::WeightedArbIndex) = wi.indexes
indextuple(wi::WeightedAdjIndex{L}) where L = ntuple(i->wi.istart+i-1, Val(L))
indextuple(wi::WeightedArbIndex{L}) where L = indexes(wi)

# Make them iterable just like numbers are
Base.iterate(x::WeightedIndex) = (x, nothing)
Base.iterate(x::WeightedIndex, ::Any) = nothing
Base.isempty(x::WeightedIndex) = false
Base.length(x::WeightedIndex) = 1

# Supporting arithmetic on the weights allows one to perform pre-scaling of gradient coefficients
Base.:(*)(wi::WeightedAdjIndex, x::Number) = WeightedAdjIndex(wi.istart, wi.weights .* x)
Base.:(/)(wi::WeightedAdjIndex, x::Number) = WeightedAdjIndex(wi.istart, wi.weights ./ x)
Base.:(*)(wi::WeightedArbIndex, x::Number) = WeightedArbIndex(wi.indexes, wi.weights .* x)
Base.:(/)(wi::WeightedArbIndex, x::Number) = WeightedArbIndex(wi.indexes, wi.weights ./ x)

### Indexing with WeightedIndex

# We inject indexing with `WeightedIndex` at a non-exported point in the dispatch heirarchy.
# This is to avoid ambiguities with methods that specialize on the array type rather than
# the index type.
Base.to_indices(A, I::Tuple{Vararg{Union{Int,WeightedIndex}}}) = I
@propagate_inbounds Base._getindex(::IndexLinear, A::AbstractVector, i::Int) = getindex(A, i)  # ambiguity resolution
@inline function Base._getindex(::IndexStyle, A::AbstractArray{T,N}, I::Vararg{Union{Int,WeightedIndex},N}) where {T,N}
    interp_getindex(A, I, ntuple(d->0, Val(N))...)
end

# The non-generated version is currently disabled due to https://github.com/JuliaLang/julia/issues/29117
# # This follows a "move processed indexes to the back" strategy, so J contains the yet-to-be-processed
# # indexes and I all the processed indexes.
# interp_getindex(A::AbstractArray{T,N}, J::Tuple{Int,Vararg{Any,L}}, I::Vararg{Int,M}) where {T,N,L,M} =
#     interp_getindex(A, Base.tail(J), I..., J[1])
# function interp_getindex(A::AbstractArray{T,N}, J::Tuple{WeightedIndex,Vararg{Any,L}}, I::Vararg{Int,M}) where {T,N,L,M}
#     wi = J[1]
#     interp_getindex1(A, indexes(wi), weights(wi), Base.tail(J), I...)
# end
# interp_getindex(A::AbstractArray{T,N}, ::Tuple{}, I::Vararg{Int,N}) where {T,N} =   # termination
#     @inbounds A[I...]  # all bounds-checks have already happened
#
# ## Handle expansion of a single dimension
# # version for WeightedAdjIndex
# @inline interp_getindex1(A, i::Int, weights::NTuple{K,Any}, rest, I::Vararg{Int,M}) where {M,K} =
#     weights[1] * interp_getindex(A, rest, I..., i) + interp_getindex1(A, i+1, Base.tail(weights), rest, I...)
# @inline interp_getindex1(A, i::Int, weights::Tuple{Any}, rest, I::Vararg{Int,M}) where M =
#     weights[1] * interp_getindex(A, rest, I..., i)
# interp_getindex1(A, i::Int, weights::Tuple{}, rest, I::Vararg{Int,M}) where M =
#     error("exhausted the weights, this should never happen")
#
# # version for WeightedArbIndex
# @inline interp_getindex1(A, indexes::NTuple{K,Int}, weights::NTuple{K,Any}, rest, I::Vararg{Int,M}) where {M,K} =
#     weights[1] * interp_getindex(A, rest, I..., indexes[1]) + interp_getindex1(A, Base.tail(indexes), Base.tail(weights), rest, I...)
# @inline interp_getindex1(A, indexes::Tuple{Int}, weights::Tuple{Any}, rest, I::Vararg{Int,M}) where M =
#     weights[1] * interp_getindex(A, rest, I..., indexes[1])
# interp_getindex1(A, indexes::Tuple{}, weights::Tuple{}, rest, I::Vararg{Int,M}) where M =
#     error("exhausted the weights and indexes, this should never happen")

@inline interp_getindex(A::AbstractArray{T,N}, J::Tuple{Int,Vararg{Any,K}}, I::Vararg{Int,N}) where {T,N,K} =
    interp_getindex(A, Base.tail(J), Base.tail(I)..., J[1])
@generated function interp_getindex(A::AbstractArray{T,N}, J::Tuple{WeightedAdjIndex{L,W},Vararg{Any,K}}, I::Vararg{Int,N}) where {T,N,K,L,W}
    ex = :(w[1]*interp_getindex(A, Jtail, Itail..., j))
    for l = 2:L
        ex = :(w[$l]*interp_getindex(A, Jtail, Itail..., j+$(l-1)) + $ex)
    end
    quote
        $(Expr(:meta, :inline))
        Jtail = Base.tail(J)
        Itail = Base.tail(I)
        j, w = J[1].istart, J[1].weights
        $ex
    end
end
@generated function interp_getindex(A::AbstractArray{T,N}, J::Tuple{WeightedArbIndex{L,W},Vararg{Any,K}}, I::Vararg{Int,N}) where {T,N,K,L,W}
    ex = :(w[1]*interp_getindex(A, Jtail, Itail..., ij[1]))
    for l = 2:L
        ex = :(w[$l]*interp_getindex(A, Jtail, Itail..., ij[$l]) + $ex)
    end
    quote
        $(Expr(:meta, :inline))
        Jtail = Base.tail(J)
        Itail = Base.tail(I)
        ij, w = J[1].indexes, J[1].weights
        $ex
    end
end
@inline interp_getindex(A::AbstractArray{T,N}, ::Tuple{}, I::Vararg{Int,N}) where {T,N} =   # termination
    @inbounds A[I...]  # all bounds-checks have already happened

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
#     itp(to_indices(itp, x)...)
# end
function gradient(itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    xi = to_indices(itp, x)
    xi == x && error("gradient of $itp not supported for position $x")
    gradient(itp, xi...)
end
@propagate_inbounds function gradient!(dest, itp::AbstractInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    dest .= gradient(itp, x...)
end
function gradient!(dest, itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    gradient!(dest, itp, to_indices(itp, x)...)
end
function hessian(itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    hessian(itp, to_indices(itp, x)...)
end
@propagate_inbounds function hessian!(dest, itp::AbstractInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    dest .= hessian(itp, x...)
end
function hessian!(dest, itp::AbstractInterpolation, x::Vararg{UnexpandedIndexTypes})
    hessian!(dest, itp, to_indices(itp, x)...)
end

# @inline function (itp::AbstractInterpolation)(x::Vararg{ExpandedIndexTypes})
#     itp.(Iterators.product(x...))
# end
# function gradient(itp::AbstractInterpolation, x::Vararg{ExpandedIndexTypes})
#     map(y->tgradient(itp, y), Iterators.product(x...))
# end
# function hessian(itp::AbstractInterpolation, x::Vararg{ExpandedIndexTypes})
#     map(y->thessian(itp, y), Iterators.product(x...))
# end
#
# tgradient(itp, y) = gradient(itp, y...)
# thessian(itp, y) = hessian(itp, y...)

# getindex is supported only for Integer indices (deliberately)
import Base: getindex
@propagate_inbounds getindex(itp::AbstractInterpolation{T,N}, i::Vararg{Integer,N}) where {T,N} = itp(i...)
@propagate_inbounds function getindex(itp::AbstractInterpolation{T,1}, i::Integer, j::Integer) where T
    @boundscheck (j == 1 || Base.throw_boundserror(itp, (i, j)))
    itp(i)
end

@inline checkbounds(::Type{Bool}, itp::AbstractInterpolation, x::Vararg{ExpandedIndexTypes,N}) where N =
    _checkbounds(BoundsCheckStyle(itp), itp, x)

_checkbounds(::CheckWillPass, itp, x) = true
_checkbounds(::NeedsCheck, itp, x) = checklubounds(lbounds(itp), ubounds(itp), x)

checklubounds(ls, us, xs) = _checklubounds(true, ls, us, xs)
_checklubounds(tf::Bool, ls, us, xs::Tuple{Number, Vararg{Any}}) =
    _checklubounds(tf & (ls[1] <= xs[1] <= us[1]), Base.tail(ls), Base.tail(us), Base.tail(xs))
_checklubounds(tf::Bool, ls, us, xs::Tuple{AbstractVector, Vararg{Any}}) =
    _checklubounds(tf & all(ls[1] .<= xs[1] .<= us[1]), Base.tail(ls), Base.tail(us), Base.tail(xs))
_checklubounds(tf::Bool, ::Tuple{}, ::Tuple{}, ::Tuple{}) = tf

maybe_clamp(itp, xs) = maybe_clamp(BoundsCheckStyle(itp), itp, xs)
maybe_clamp(::NeedsCheck, itp, xs) = clamp.(xs, lbounds(itp), ubounds(itp))
maybe_clamp(::CheckWillPass, itp, xs) = xs

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
include("gridded/gridded.jl")
include("monotonic/monotonic.jl")
include("extrapolation/extrapolation.jl")
include("scaling/scaling.jl")
include("utils.jl")
include("io.jl")
include("convenience-constructors.jl")
include("deprecations.jl")

end # module
