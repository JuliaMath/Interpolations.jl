# Similar to ExtrapDimSpec but for only a single dimension
const ExtrapSpec = Union{BoundaryCondition,Tuple{BoundaryCondition,BoundaryCondition}}

# Macro to get create ExtrapSpec for checking if a KnotIterator has a given BC
# for forward iteration
macro FwdExtrapSpec(bc)
    :( Union{$bc,Tuple{BoundaryCondition,$bc}} )
end
struct KnotIterator{T,ET <: ExtrapSpec}
    knots::Vector{T}
    bc::ET
    nknots::Int
    KnotIterator{T,ET}(k::AbstractArray{T}, bc::ET) where {T,ET} = new(k, bc, length(k))
end

# Dispatch on outputs of k = getknots and bc = etpflag to create one
# KnotIterator per interpolation axis
function KnotIterator(k::NTuple{N,AbstractArray}, bc::NTuple{N, ExtrapSpec}) where {N}
    # One ExtrapSpec per Axis
    map(KnotIterator, k, bc)
end
function KnotIterator(k::NTuple{N,AbstractArray}, bc::BoundaryCondition) where {N}
    # All Axes use the same BoundaryCondition
    map(x -> KnotIterator(x, bc), k)
end
function KnotIterator(k::AbstractArray{T}, bc::ET) where {T,ET <: ExtrapDimSpec}
    # Prior dispatches will end up here: Axis knots: k and boundary condition bc
    KnotIterator{T,ET}(k, bc)
end

const RepeatKnots = Union{Periodic,Reflect}
IteratorSize(::Type{KnotIterator}) = SizeUnknown()
IteratorSize(::Type{KnotIterator{T}}) where {T} = SizeUnknown()

IteratorSize(::Type{KnotIterator{T,ET}}) where {T,ET} = _knot_iter_size(ET)
_knot_iter_size(::Type{<:BoundaryCondition}) = HasLength()
_knot_iter_size(::Type{<:RepeatKnots}) = IsInfinite()
_knot_iter_size(::Type{Tuple{RevBC, FwdBC}}) where {RevBC,FwdBC} = _knot_iter_size(FwdBC)

length(iter::KnotIterator) = _knot_length(iter, IteratorSize(iter))
_knot_length(iter::KnotIterator, ::HasLength) = iter.nknots
size(iter::KnotIterator) = (length(iter),)

IteratorEltype(::Type{KnotIterator{T,ET}}) where {T,ET} = HasEltype()
eltype(::Type{KnotIterator{T,ET}}) where {T,ET} = T

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
    # Constuct seperate KnotIterator for each dimension, and combine them using
    # Iterators.product
    k = getknots(itp)
    bc = Throw()
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? iter[1] : Iterators.product(iter...)
end

function knots(etp::AbstractExtrapolation)
    # Constuct seperate KnotIterator for each dimension, and combine them using
    # Iterators.product
    k = getknots(etp)
    bc = etpflag(etp)
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? iter[1] : Iterators.product(iter...)
end

# For non-repeating ET's iterate through once
iterate(iter::KnotIterator) where {T} = iterate(iter, 1)
iterate(iter::KnotIterator, idx::Integer) = idx <= iter.nknots ? (iter.knots[idx], idx+1) : nothing

# For repeating knots state is the knot index + offset value
function iterate(iter::KnotIterator{T,ET}) where {T,ET <: @FwdExtrapSpec RepeatKnots}
    iterate(iter, (1, zero(T)))
end

# Periodic: Iterate over knots, updating the offset each cycle
function iterate(iter::KnotIterator{T,ET}, state) where {T, ET <: @FwdExtrapSpec(Periodic)}
    state === nothing && return nothing
    curidx, offset = state[1], state[2]

    # Increment offset after cycling over all the knots
    if mod(curidx, iter.nknots-1) != 0
        nextstate = (curidx+1, offset)
    else
        knotrange = iter.knots[end] - iter.knots[1]
        nextstate = (curidx+1, offset+knotrange)
    end

    # Get the current knot
    knot = iter.knots[periodic(curidx, 1, iter.nknots)] + offset
    return knot, nextstate
end

# Reflect: Iterate over knots, updating the offset after a forward and backwards
# cycle
function iterate(iter::KnotIterator{T, ET}, state) where {T, ET <: @FwdExtrapSpec(Reflect)}
    state === nothing && return nothing
    curidx, offset = state[1], state[2]

    # Increment offset after a forward + backwards pass over the knots
    cycleidx = mod(curidx, 2*iter.nknots-1)
    if  cycleidx != 0
        nextstate = (curidx+1, offset)
    else
        knotrange = iter.knots[end] - iter.knots[1]
        nextstate = (curidx+1, offset+2*knotrange)
    end

    # Map global index onto knots, and get that knot
    idx = reflect(curidx, 1, iter.nknots)
    knot = iter.knots[idx]

    # Add offset to map local knot to global position
    if 0 < cycleidx <= iter.nknots
        # Forward pass
        knot = knot + offset
    else
        # backwards pass
        knot = offset + 2*iter.knots[end] - knot
    end

    return knot, nextstate
end
