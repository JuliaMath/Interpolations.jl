# Similar to ExtrapDimSpec but for only a single dimension
const ExtrapSpec = Union{BoundaryCondition,Tuple{BoundaryCondition,BoundaryCondition}}

# Type Alias to get Boundary Condition or forward boundary conditions if
# directional
const FwdExtrapSpec{FwdBC} = Union{FwdBC, Tuple{BoundaryCondition, FwdBC}}
const RevExtrapSpec{RevBC} = Union{RevBC, Tuple{RevBC, BoundaryCondition}}

"""
    KnotIterator{T,ET}(k::AbstractArray{T}, bc::ET)

Defines an iterator over the knots in `k` based on the boundary conditions `bc`.

# Fields
- `knots::Vector{T}` The interpolated knots of the axis to iterate over
- `bc::ET` The Boundary Condition for the axis
- `nknots::Int` The number of interpolated knots (ie. `length(knots)`)

`ET` is `Union{BoundaryCondition,Tuple{BoundaryCondition,BoundaryCondition}}`

# Iterator Interface
The following methods defining Julia's iterator interface have been defined

`IteratorSize(::Type{KnotIterator})` -> Will return one of the following
- `Base.IsInfinite` if the iteration will produces an infinite sequence of knots
- `Base.HasLength` if iteration will produce a finite sequence of knots
- `Base.SizeUnknown` if we can't decided from only the type information

`length` and `size` -> Are defined if IteratorSize is HasLength, otherwise will
raise a MethodError.

`IteratorEltype` will always return `HasEltype`, as we always track the data
types of the knots

`eltype` will return the data type of the knots

`iterate` Defines iteration over the knots starting from the first one and
moving in the forward direction along the axis.

# Knots for Multi-dimensional Interpolants
Iteration over the knots of a multi-dimensional interpolant is done by wrapping
multiple KnotIterator within `Iterators.product`.

"""
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

Iterator will yield scalar values for interpolations over a single dimension,
and tuples of coordinates for higher dimension interpolations. Iteration over
higher dimensions is taken as the product of knots along each dimension.

ie. Iterator.product(knots on first dim, knots on 2nd dim,...)

Extrapolations with Periodic or Reflect boundary conditions, will produce an
infinite sequence of knots.

# Example
```jldoctest
julia> using Interpolations;

julia> etp = LinearInterpolation([1.0, 1.2, 2.3, 3.0], rand(4); extrapolation_bc=Periodic());

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
    # Construct separate KnotIterator for each dimension, and combine them using
    # Iterators.product
    k = getknots(itp)
    bc = Throw()
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? iter[1] : Iterators.product(iter...)
end

function knots(etp::AbstractExtrapolation)
    # Construct separate KnotIterator for each dimension, and combine them using
    # Iterators.product
    k = getknots(etp)
    bc = etpflag(etp)
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? iter[1] : Iterators.product(iter...)
end

# Initialize iteration based on KnotIterator
_knot_start(::KnotIterator) = 1
_knot_start(::KnotIterator{T,ET}) where {T,ET <: FwdExtrapSpec{RepeatKnots}} = (1, zero(T))

# For non-repeating ET's iterate through once
iterate(iter::KnotIterator) = iterate(iter, _knot_start(iter))
iterate(iter::KnotIterator, idx::Integer) = idx <= iter.nknots ? (iter.knots[idx], idx+1) : nothing

# Periodic: Iterate over knots, updating the offset each cycle
function iterate(iter::KnotIterator{T,ET}, state::Tuple) where {T, ET <: FwdExtrapSpec{Periodic}}
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
function iterate(iter::KnotIterator{T, ET}, state) where {T, ET <: FwdExtrapSpec{Reflect}}
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
struct KnotRange{T, L <: IteratorSize}
    iter::KnotIterator{T}
    start
    stop
    function KnotRange(iter::KnotIterator{T}, start, stop) where {T,ET}
        L = IteratorSize(iter) == HasLength() || stop !== nothing ? HasLength : IsInfinite
        new{T, L}(iter, start, stop)
    end
end

#KnotRange(iter::KnotIterator{T}, start, stop) where {T} = KnotRange{T}(iter, start, stop)

IteratorSize(::Type{KnotRange}) = SizeUnknown()
IteratorSize(::Type{KnotRange{T}}) where {T} = SizeUnknown()
IteratorSize(::Type{<:KnotRange{T,L}}) where {T,L} = L()
IteratorEltype(::Type{<:KnotRange}) = HasEltype()
eltype(::Type{<:KnotRange{T}}) where {T} = T

function length(iter::KnotRange{T,HasLength}) where {T}
    sdx = _knot_start(iter.iter, iter.start)
    edx = _knot_stop(iter.iter, iter.stop)
    _knot_iter_length(sdx, edx)
end
_knot_iter_length(sdx::Int, edx::Int) = max(edx - sdx + 1, 0)
_knot_iter_length(sdx, edx) = _knot_iter_length(first(sdx), first(edx))
_knot_iter_length(::Any, ::Nothing) = 0
_knot_iter_length(::Nothing, ::Any) = 0

size(iter::KnotRange{T,HasLength}) where {T} = (length(iter),)

knotsbetween(itp; start=nothing, stop=nothing) = knotsbetween(itp, start, stop)
knotsbetween(iter::KnotIterator, start, stop) = KnotRange(iter, start, stop)

# If AbstractInterpolation -> Get Knots -> Then Convert to Knot Range
knotsbetween(itp::AbstractInterpolation, start, stop) = knotsbetween(knots(itp), start, stop)

knotsbetween(::KnotIterator, ::Nothing, ::Nothing) =
    throw(ArgumentError("At least one of `start` or `stop` must be specified"))

# Expand Nothing on Start/Stop if tuples
function knotsbetween(iter::Iterators.ProductIterator, ::Nothing, stop::Tuple)
    N = length(stop)
    knotsbetween(iter, Tuple(repeat([nothing], N)), stop)
end
function knotsbetween(iter::Iterators.ProductIterator, start::Tuple, ::Nothing)
    N = length(start)
    knotsbetween(iter, start, Tuple(repeat([nothing], N)))
end

function knotsbetween(iter::Iterators.ProductIterator, start::Tuple, stop::Tuple)
    # Get KnotIterator for each dimension -> Wrap with KnotRange -> Recombine
    # with Iterators.product
    kiter = map(knotsbetween, iter.iterators, start, stop)
    Iterators.product(kiter...)
end

function iterate(iter::KnotRange, state)
    y = iterate(iter.iter, state)
    y === nothing && return nothing
    if iter.stop !== nothing
        return y[1] < iter.stop ? y : nothing
    else
        return y
    end
end

# If state isn't provided -> Compute it for the KnotRange
iterate(iter::KnotRange) = iterate(iter, _knot_start(iter.iter, iter.start))

# If start is nothing -> Fall back to default start
_knot_start(iter::KnotRange, start) = _knot_start(iter.iter, start)
_knot_start(iter::KnotIterator, ::Nothing) = _knot_start(iter)

_knot_stop(iter::KnotRange, stop) = _knot_stop(iter.iter, stop)
_knot_stop(iter::KnotIterator, ::Nothing) = _knot_stop(iter)
_knot_stop(iter::KnotIterator) = length(iter)

etptypes(ET::Type{<:BoundaryCondition}) = (ET, ET)
etptypes(::Type{Tuple{RevBC, FwdBC}}) where {RevBC,FwdBC} = (RevBC, FwdBC)

# If start is provided, decided if interpolated knot, or left / right
# extrapolated knot
function _knot_start(iter::KnotIterator{T,ET}, start) where {T, ET <: ExtrapDimSpec}
    if iter.knots[1] <= start
        # Starting knot is interpolated / On the Right Side
        return _knot_start(etptypes(ET)[2], iter, start)
    else
        # Starting knot is extrapolated on the left side
        return _knot_start(etptypes(ET)[1], iter, start)
    end
end
function _knot_stop(iter::KnotIterator{T,ET}, stop) where {T, ET <: ExtrapDimSpec}
    if stop <= iter.knots[end]
        # Ending knot is interpolated / On the left Side
        return _knot_stop(etptypes(ET)[1], iter, stop)
    else
        # Starting knot is extrapolated on the right side
        return _knot_stop(etptypes(ET)[2], iter, stop)
    end
end

# Starting iterator state for a non-repeating knot
function _knot_start(::Type{<:BoundaryCondition}, iter, start)
    findfirst(start .< iter.knots)
end
function _knot_stop(::Type{<:BoundaryCondition}, iter, stop)
    findlast(iter.knots .< stop)
end

# Starting iterator state for a periodic knots
function _knot_start(::Type{<:Periodic}, iter::KnotIterator{T}, start) where {T}
    # Find starting offset
    knotrange = iter.knots[end] - iter.knots[1]
    cycle = floor(Int, (start-iter.knots[1])/ knotrange)

    # Find starting index
    cyclepos = periodic(start, iter.knots[1], iter.knots[end])
    cycleidx = findfirst(cyclepos .< iter.knots)::Int
    if cycleidx === iter.nknots
        cycleidx = 1
        cycle += 1
    end

    offset = knotrange * cycle
    idx = (iter.nknots-1) * cycle + cycleidx

    return idx, offset
end
function _knot_stop(::Type{<:Periodic}, iter::KnotIterator{T}, stop) where {T}
    # Find stopping offset
    knotrange = iter.knots[end] - iter.knots[1]
    cycle = floor(Int, (stop-iter.knots[1])/ knotrange)

    # Find stopping index
    cyclepos = periodic(stop, iter.knots[1], iter.knots[end])
    if cyclepos == iter.knots[1]
        cycleidx = iter.nknots -1
        cycle -= 1
    else
        cycleidx = findlast(iter.knots .< cyclepos)::Int
    end
    offset = knotrange * cycle
    idx = (iter.nknots-1) * cycle + cycleidx

    return idx, offset
end

# Starting iterator state for a reflecting knots
function _knot_start(::Type{<:Reflect}, iter::KnotIterator{T}, start) where {T}
    # Find starting offset
    knotrange = iter.knots[end] - iter.knots[1]
    cycleidx = floor(Int, (start-iter.knots[1])/ knotrange)
    offset = knotrange * cycleidx

    # Find starting index
    inbound_start = reflect(start, iter.knots[1], iter.knots[end])
    idx = iter.nknots * cycleidx + findfirst(inbound_start .< iter.knots)

    return idx, offset
end
function _knot_stop(::Type{<:Reflect}, iter::KnotIterator{T}, stop) where {T}
    # Find stopping offset
    knotrange = iter.knots[end] - iter.knots[1]
    cycleidx = floor(Int, (stop-iter.knots[1])/ knotrange)
    offset = knotrange * cycleidx

    # Find stopping index
    inbound_stop = reflect(stop, iter.knots[1], iter.knots[end])
    idx = iter.nknots * cycleidx + findlast(iter.knots .< inbound_stop)

    return idx, offset
end
