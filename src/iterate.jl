
export knots

# Similar to ExtrapDimSpec but for only a single dimension
const ExtrapSpec = Union{BoundaryCondition,Tuple{BoundaryCondition,BoundaryCondition}}

struct KnotIterator{T,ET <: ExtrapSpec}
    knots::Vector{T}
    bc::ET
    nknots::Int
    KnotIterator{T,ET}(k::AbstractArray{T}, bc::ET) where {T,ET} = new(k, bc, length(k))
end

KnotIterator(k::Tuple, bc::Tuple) = map(KnotIterator, k, bc)
KnotIterator(k::Tuple, bc::BoundaryCondition) = map(x -> KnotIterator(x, bc), k)
KnotIterator(k::AbstractArray{T}, bc::ET) where {T,ET <: ExtrapSpec} = KnotIterator{T,ET}(k, bc)

const RepeatKnots = Union{Periodic,Reflect}
Base.length(iter::KnotIterator{T,ET}) where {T,ET} = length(iter.knots)
Base.IteratorSize(::Type{KnotIterator{T,ET}}) where {T,ET <: RepeatKnots} = Base.IsInfinite()
Base.length(::KnotIterator{T,ET}) where {T,ET <: RepeatKnots} = Int(Inf)

Base.IteratorEltype(::Type{KnotIterator{T,ET}}) where {T,ET} = Base.HasEltype()
Base.eltype(::Type{KnotIterator{T,ET}}) where {T,ET} = T

"""
    knots(itp::AbstractInterpolation)
    knots(etp::AbstractExtrapolation)

Returns an iterator over knot locations for an AbstractInterpolation or
AbstractExtrapolation.

Iterator will yeild scalar values for interpolations over a single dimension,
and tuples of coordinates for higher dimension interpolations. Iteration over
higher dimension is taken as the product of knots on each dimensions.

ie. Iterator.product(knots on first dim, knots on 2nd dim,...)

Extrapolations with Periodic or Reflect boundary conditions, will produce an
infinite sequence of knots.

# Example
```julia-repl
julia> etp = LinearInterpolation([1.0, 1.2, 2.3, 3.0], rand(4); extrapolation_bc=Periodic())
julia> Iterators.take(knots(etp), 5) |> collect
5-element Array{Float64,1}:
 1.0
 1.2
 2.3
 3.0
 3.2
```
"""
function knots(itp::AbstractInterpolation)
    # For interpolations -> Return an array or Iterators.flatten of the knots
    I = Iterators.product(getknots(itp)...)
    ndims(itp) == 1 ? Iterators.flatten(I) : I
end

function knots(etp::AbstractExtrapolation)
    # Constuct seperate KnotIterator for each dimension, and combine them using
    # Iterators.product
    k = getknots(etp)
    bc = etpflag(etp)
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? only(iter) : Iterators.ProductIterator(iter)
end

# Start at the first knot, with zero offset by default
Base.iterate(iter::KnotIterator{T}) where {T} = iterate(iter, (1, zero(T)))

# Iterate over knots until curidx is nothing, dispatching to nextstate to handle
# repeated / non-repeating, and directional boundary conditions
function Base.iterate(iter::KnotIterator, state)
    curidx, offset = state[1], state[2]
    isnothing(curidx) && return nothing
    idx = inbounds_index(iter.bc, (1, iter.nknots), curidx, iter, state)
    iter.knots[idx] + offset, nextstate(iter, state)
end

# For non-repeating knots, iterate over knots then return nothing
function nextstate(iter::KnotIterator, state)
    idx, offset = state
    nextidx = idx < iter.nknots ? idx + 1 : nothing
    (nextidx, offset)
end

# For repeating knots iterate over knots, then increment offset
# The last knot is "skipped" as it has the same coordinate as the first knot
function nextstate(iter::KnotIterator{T,ET}, state) where {T,ET<:RepeatKnots}
    idx, offset = state
    if idx + 1 < iter.nknots
        return idx + 1, offset
    else
        knotrange = iter.knots[end] - iter.knots[1]
        return 1, offset + knotrange
    end
end
