module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,
    scale,
    bounds,

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

    constant_interpolation,
    linear_interpolation,
    cubic_spline_interpolation,

    knots,
    knotsbetween

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl
    # monotonic/monotonic.jl
    # scaling/scaling.jl
    # hermite/cubic.jl

using LinearAlgebra, SparseArrays
using StaticArrays, WoodburyMatrices, Ratios, AxisAlgorithms, OffsetArrays
using ChainRulesCore

using Base: @propagate_inbounds, HasEltype, EltypeUnknown, HasLength, IsInfinite,
    SizeUnknown, Indices
import Base: convert, size, axes, promote_rule, ndims, eltype, checkbounds, axes1,
    iterate, length, IteratorEltype, IteratorSize, firstindex, getindex, LogicalIndex

abstract type Flag end
abstract type InterpolationType <: Flag end
"`NoInterp()` indicates that the corresponding axis must use integer indexing (no interpolation is to be performed)"
struct NoInterp <: InterpolationType end
abstract type GridType <: Flag end
"`OnGrid()` indicates that the boundary condition applies at the first & last nodes"
struct OnGrid <: GridType end
"`OnCell()` indicates that the boundary condition applies a half-gridspacing beyond the first & last nodes"
struct OnCell <: GridType end

const DimSpec{T} = Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract type AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType}} <: AbstractArray{T,N} end
abstract type AbstractInterpolationWrapper{T,N,ITPT,IT} <: AbstractInterpolation{T,N,IT} end
abstract type AbstractExtrapolation{T,N,ITPT,IT} <: AbstractInterpolationWrapper{T,N,ITPT,IT} end

function Base.:(==)(itp1::AbstractInterpolation, itp2::AbstractInterpolation)
    propertynames(itp1) == propertynames(itp2) || return false
    for pname in propertynames(itp1)
        getproperty(itp1, pname) == getproperty(itp2, pname) || return false
    end
    return true
end

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

where `gt` is the grid type, e.g., [`OnGrid()`](@ref) or [`OnCell()`](@ref).
"""
abstract type BoundaryCondition <: Flag end
# Put the gridtype into the boundary condition, since that's all it affects (see issue #228)
# Nothing is used for extrapolation
"`Throw(gt)` causes beyond-the-edge extrapolation to throw an error"
struct Throw{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`Flat(gt)` sets the extrapolation slope to zero"
struct Flat{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`Line(gt)` uses a constant slope for extrapolation"
struct Line{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`Free(gt)` the free boundary condition makes sure the interpoland has a continuous third derivative at the second-to-outermost cell boundary"
struct Free{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`Periodic(gt)` applies periodic boundary conditions"
struct Periodic{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`Reflect(gt)` applies reflective boundary conditions"
struct Reflect{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"`InPlace(gt)` is a boundary condition that allows prefiltering to occur in-place (it typically requires padding)"
struct InPlace{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
"""`InPlaceQ(gt)` is similar to `InPlace(gt)`, but is exact when the values being interpolated
arise from an underlying quadratic. It is primarily useful for testing purposes,
allowing near-exact (to machine precision) comparisons against ground truth."""
struct InPlaceQ{GT<:Union{GridType,Nothing}} <: BoundaryCondition gt::GT end
const Natural = Line

(::Type{BC})() where BC<:BoundaryCondition = BC(nothing)
for BC in (:Throw, :Flat, :Line, :Free, :Reflect, :InPlace, :InPlaceQ, :Periodic)
    eval(quote
        $BC(::Type{GT}) where GT <: GridType = $BC{GT}(GT())
        $BC{GT}(::Type{GT}) where GT <: GridType = $BC{GT}(GT())
    end)
end
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

"""
    bounds(itp::AbstractInterpolation)

Return the `bounds` of the domain of `itp` as a tuple of `(min, max)` pairs for each coordinate. This is best explained by example:

```jldoctest
julia> itp = interpolate([1 2 3; 4 5 6], BSpline(Linear()));

julia> bounds(itp)
((1, 2), (1, 3))

julia> data = 1:3;

julia> knots = ([10, 11, 13.5],);

julia> itp = interpolate(knots, data, Gridded(Linear()));

julia> bounds(itp)
((10.0, 13.5),)
```
"""
bounds(itp::AbstractInterpolation) = map(twotuple, lbounds(itp), ubounds(itp))
bounds(itp::AbstractInterpolation, d) = bounds(itp)[d]

itptype(::Type{<:AbstractInterpolation{T,N,IT}}) where {T,N,IT<:DimSpec{InterpolationType}} = IT
itptype(itp::AbstractInterpolation) = itptype(typeof(itp))
# The following have been defined by AbstractArray
# ndims(::Type{<:AbstractInterpolation{T,N,<:DimSpec{InterpolationType}}}) where {T,N} = N
# ndims(itp::AbstractInterpolation) = ndims(typeof(itp))
# eltype(::Type{<:AbstractInterpolation{T,N,<:DimSpec{InterpolationType}}}) where {T,N} = T
# eltype(itp::AbstractInterpolation) = eltype(typeof(itp))

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
```julia
ret = (1-f)*a[i] + f*a[i+1]
```
where `f = x-i` lies between 0 and 1. This can be represented as `a[wi]`, where
```julia
wi = WeightedIndex(i:i+1, (1-f, f))
```
i.e.,
```julia
ret = sum(a[indexes] .* weights)
```

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
```julia
gwi = WeightedIndex(i:i+1, (-1, 1))
g1 = a[gwi]
```
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
    @noinline mismatch(indexes, weights) =
        throw(ArgumentError("the length of indexes must match weights, got $indexes vs $weights"))
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

# We inject `WeightedIndex` as a non-exported indexing point with a `InterpGetindex` wrapper.
# `InterpGetindex` is not a subtype of `AbstractArray`. This ensures that the overload applies to all array types.
struct InterpGetindex{N,A<:AbstractArray{<:Any,N}}
    coeffs::A
    InterpGetindex(itp::AbstractInterpolation) = InterpGetindex(coefficients(itp))
    InterpGetindex(A::AbstractArray) = new{ndims(A),typeof(A)}(A)
end
@inline Base.getindex(A::InterpGetindex{N}, I::Vararg{Union{Int,WeightedIndex},N}) where {N} =
    interp_getindex(A.coeffs, ntuple(zero, Val(N)), map(indexflag, I)...)
@inline indexflag(I) = indextuple(I), weights(I)

# Direct recursion would allow more eager inference before julia 1.11.
# Normalize all index into the same format.
struct One end  # Singleton for express weights of no-interp dims
indextuple(I::Int) = (I,)
weights(::Int) = (One(),)

struct Zero end # Singleton for dim expansion termination

# A recursion-based `interp_getindex`, which follows a "move processed indexes to the back" strategy
# `I` contains the processed index, and (wi1, wis...) contains the yet-to-be-processed indexes
# Here we handle the expansion of a single dimension.
@inline function interp_getindex(A, I, (is, ws)::NTuple{2,Tuple}, wis...)
    itped1 = interp_getindex(A, (Base.tail(I)..., is[end]), wis...)
    witped = interp_getindex(A, I, (Base.front(is), Base.front(ws)), wis...)
    _weight_itp(ws[end], itped1, witped)
end
interp_getindex(_, _, ::NTuple{2,Tuple{}}, ::Vararg) = Zero()
# Termination
@inline interp_getindex(A::AbstractArray{T,N}, I::Dims{N}) where {T,N} =
    @inbounds A[I...] # all bounds-checks have already happened

_weight_itp(w, i, wir) = w * i + wir
_weight_itp(::One, i, ::Zero) = i
_weight_itp(w, i, ::Zero) = w * i

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
# Explicitly handle boolean vectors to avoid ambiguity with Base.checkbounds
@inline checkbounds(::Type{Bool}, itp::AbstractInterpolation, x::AbstractVector{Bool}) =
    _checkbounds(BoundsCheckStyle(itp), itp, x)

@inline checkbounds(::Type{Bool}, itp::AbstractInterpolation, x::LogicalIndex) =
    _checkbounds(BoundsCheckStyle(itp), itp, x)

@inline checkbounds(::Type{Bool}, itp::AbstractInterpolation, x::LogicalIndex{<:Any,<:AbstractArray{Bool,1}}) =
    _checkbounds(BoundsCheckStyle(itp), itp, x)

_checkbounds(::CheckWillPass, itp, x) = true
_checkbounds(::NeedsCheck, itp, x) = checklubounds(lbounds(itp), ubounds(itp), x)

checklubounds(ls, us, xs) = _checklubounds(true, ls, us, xs)
_checklubounds(tf::Bool, ls::Tuple, us::Tuple, xs::Tuple) =
    _checklubounds(tf & allbetween(ls[1], xs[1], us[1]), Base.tail(ls), Base.tail(us), Base.tail(xs))
_checklubounds(tf::Bool, ::Tuple{}, ::Tuple{}, xs::Tuple) =
    _checklubounds(tf & all(isone, xs[1]), (), (), Base.tail(xs))
_checklubounds(tf::Bool, ls::Tuple, us::Tuple, ::Tuple{}) =
    _checklubounds(tf & (ls[1] == us[1]), Base.tail(ls), Base.tail(us), ())
_checklubounds(tf::Bool, ::Tuple{}, ::Tuple{}, ::Tuple{}) = tf

maybe_clamp(itp, xs) = maybe_clamp(BoundsCheckStyle(itp), itp, xs)
maybe_clamp(::NeedsCheck, itp, xs) = map(clamp, xs, lbounds(itp), ubounds(itp))
maybe_clamp(::CheckWillPass, itp, xs) = xs

# this strips arbitrary layers of ForwardDiff.Dual, returning the innermost value
# it's other methods are defined in InterpolationsForwardDiffExt.jl
just_dual_value(x::Number) = x

Base.hash(x::AbstractInterpolation, h::UInt) = Base.hash_uint(3h - objectid(x))
Base.hash(x::AbstractExtrapolation, h::UInt) = Base.hash_uint(3h - objectid(x))

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
include("lanczos/lanczos.jl")
include("lanczos/lanczos_opencv.jl")
include("iterate.jl")
include("chainrules/chainrules.jl")
include("hermite/cubic.jl")
include("gpu_support.jl")

end # module
