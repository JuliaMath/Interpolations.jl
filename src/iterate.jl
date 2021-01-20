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

# For non-repeating ET's iterate through once
firstindex(::KnotIterator) = 1
iterate(iter::KnotIterator) = iterate(iter, firstindex(iter))
function iterate(iter::KnotIterator{T,ET}, idx) where {T, ET}
    if checkbounds(Bool, iter, idx)
        iter[idx], idx+1
    else
        nothing
    end
end

splitExtrapDimSpec(::Type{ET}) where {ET <: BoundaryCondition} = (ET, ET)
splitExtrapDimSpec(::Type{<:Tuple{L,R}}) where {L,R} = (L,R)

function checkbounds(::Type{Bool}, iter::KnotIterator{T,ET}, idx::Int) where {T,ET<:ExtrapDimSpec}
    # Get Left/Right Extrapolation Boundary Conditions for the KnotIterator
    left, right = splitExtrapDimSpec(ET)

    # Check Left/Right Boundary Limits
    leftcheck = left <: RepeatKnots ? true : 1 <= idx
    rightcheck = right <: RepeatKnots ? true : idx <= iter.nknots
    leftcheck && rightcheck
end

# Raise BoundsError if knots are not extrapolated
function getindex(iter::KnotIterator{T, ET}, idx::Int) where {T, ET <: ExtrapDimSpec}
    # Get Left/Right Extrapolation Boundary Conditions for the KnotIterator
    left, right = splitExtrapDimSpec(ET)

    # Construct function to call the correct getknotindex if idx is to the left,
    # inside, or right of the interpolated knots
    if idx < 1
        getknotindex(left, iter, idx)::T
    elseif idx <= iter.nknots
        iter.knots[idx]
    else
        getknotindex(right, iter, idx)::T
    end
end

# For Non-repeating Knots -> Raise BoundsError
function getknotindex(::Type{<:BoundaryCondition}, iter, idx)
    throw(BoundsError(iter, idx))
end

# Get Knots using a Periodic Boundary Condition
@inline function getknotindex(::Type{<:Periodic}, iter, idx)
    knotrange = iter.knots[end] - iter.knots[1]
    offset = knotrange * (cld(idx, iter.nknots-1)-1)
    iter.knots[periodic(idx, 1, iter.nknots)] + offset
end

# Get Knots using a Reflect Boundary Condition
@inline function getknotindex(::Type{<:Reflect}, iter, idx)
    # Get the current knot
    knotrange = iter.knots[end] - iter.knots[1]
    offset = 2*knotrange * (cld(idx, 2*iter.nknots-2)-1)
    cycleidx = mod(idx-1, 2*iter.nknots-2)+1

    if iter.nknots < cycleidx
        offset + 2*iter.knots[end] - iter.knots[2*iter.nknots-cycleidx]
    else
        iter.knots[cycleidx] + offset
    end
end

function nextknotidx(iter::KnotIterator{T,ET}, x) where {T, ET}
    left, right = splitExtrapDimSpec(ET)
    if x < iter.knots[end]
        # Next knot is inbounds or to the left
        nextknotidx(left, iter.knots, x)
    else
        # Next knot is to the right
        nextknotidx(right, iter.knots, x)
    end
end

function priorknotidx(iter::KnotIterator{T,ET}, x) where {T, ET}
    left, right = splitExtrapDimSpec(ET)
    if x <= iter.knots[1]
        # Prior knot is to the left
        priorknotidx(left, iter.knots, x)
    else
        # prior knot is inbounds or to the right
        priorknotidx(right, iter.knots, x)
    end
end

function nextknotidx(::Type{<:BoundaryCondition}, knots::Vector, start)
    knots[end] <= start ? length(knots)+1 : findfirst(start .< knots)
end
function priorknotidx(::Type{<:BoundaryCondition}, knots::Vector, stop)
    stop <= knots[1] ? 0 : findlast(knots .< stop)
end

# Starting iterator state for a periodic knots
function nextknotidx(::Type{<:Periodic}, knots::Vector, start)
    # Find starting offset
    knotrange = knots[end] - knots[1]
    cycle = floor(Int, (start-knots[1])/ knotrange)

    # Find starting index
    cyclepos = periodic(start, knots[1], knots[end])
    cycleidx = findfirst(cyclepos .< knots)::Int
    if cycleidx === length(knots)
        cycleidx = 1
        cycle += 1
    end

    idx = (length(knots)-1) * cycle + cycleidx

    return idx
end
function priorknotidx(::Type{<:Periodic}, knots::Vector, stop)
    # Find stopping offset
    knotrange = knots[end] - knots[1]
    cycle = floor(Int, (stop-knots[1])/ knotrange)

    # Find stopping index
    cyclepos = periodic(stop, knots[1], knots[end])
    if cyclepos == knots[1]
        cycleidx = length(knots) -1
        cycle -= 1
    else
        cycleidx = findlast(knots .< cyclepos)::Int
    end
    idx = (length(knots)-1) * cycle + cycleidx

    return idx
end

# Starting iterator state for a reflecting knots
function nextknotidx(::Type{<:Reflect}, knots::Vector, start)
    refknots = knots[end] .+ (cumsum∘reverse∘diff)(knots)
    knots = vcat(knots, refknots)
    nextknotidx(Periodic, knots, start)
end
function priorknotidx(::Type{<:Reflect}, knots::Vector, stop)
    refknots = knots[end] .+ (cumsum∘reverse∘diff)(knots)
    knots = vcat(knots, refknots)
    priorknotidx(Periodic, knots, stop)
end

struct KnotRange{T, R}
    iter::KnotIterator{T}
    range::R
    KnotRange(iter::KnotIterator{T}, range::R) where {T,R} = new{T,R}(iter, range)
end

function KnotRange(iter::KnotIterator, start, stop)
    # Get the index of the first knot larger than start or firstindex(iter) iff
    # start is missing
    sdx::Int = start === nothing ? firstindex(iter) : nextknotidx(iter, start)

    # Generate knot index range from provided stop
    if stop === nothing
        if IteratorSize(iter) === IsInfinite()
            # Stop is missing and iter generates infinite knots
            range = Iterators.countfrom(sdx)
        else
            # Stop is missing and iter generates a finite number of knots
            range = sdx:iter.nknots
        end
    else
        # Stop is provided -> Iterate from start to stop
        edx = priorknotidx(iter, stop)::Int
        range = sdx:edx
    end
    KnotRange(iter, range)
end

# Iterator Size is Unknown until we knot range's type
IteratorSize(::Type{KnotRange}) = SizeUnknown()
IteratorSize(::Type{KnotRange{T}}) where {T} = SizeUnknown()
IteratorSize(::Type{KnotRange{T,R}}) where {T, R <: UnitRange} = HasLength()
IteratorSize(::Type{KnotRange{T,R}}) where {T, R <: Iterators.Count} = IsInfinite()

IteratorEltype(::Type{<:KnotRange}) = HasEltype()
eltype(::Type{<:KnotRange{T}}) where {T} = T

# Dispatch length and size to range
length(iter::KnotRange) = length(iter.range)
size(iter::KnotRange) = size(iter.range)

knotsbetween(itp; start=nothing, stop=nothing) = knotsbetween(itp, start, stop)

# If AbstractInterpolation -> Get Knots -> Then Convert to Knot Range
knotsbetween(itp::AbstractInterpolation, start, stop) = knotsbetween(knots(itp), start, stop)

knotsbetween(::KnotIterator, ::Nothing, ::Nothing) =
    throw(ArgumentError("At least one of `start` or `stop` must be specified"))

knotsbetween(iter::KnotIterator, start, stop) = KnotRange(iter, start, stop)

# Expand Nothing on Start/Stop if tuples
function knotsbetween(iter::Iterators.ProductIterator, start::Union{Nothing, Tuple}, stop::Union{Nothing, Tuple})
    kiter = map_ranges(iter.iterators, start, stop)
    Iterators.product(kiter...)
end

# Map iter, start and stop to N KnotRanges, handling start/stop being nothing
function map_ranges(iter, start, ::Nothing)
    map((i, s) -> knotsbetween(i, s, nothing), iter, start)
end
function map_ranges(iter, ::Nothing, stop)
    map((i, s) -> knotsbetween(i, nothing, s), iter, stop)
end
function map_ranges(iter, start, stop)
    map(knotsbetween, iter, start, stop)
end

iterate(iter::KnotRange) = iterate(iter, first(iter.range))

function iterate(iter::KnotRange{T,R}, state) where {T, R <: UnitRange}
    last(iter.range) < state && return nothing
    iter.iter[state], state+1
end

function iterate(iter::KnotRange{T,R}, state) where {T, R <: Iterators.Count}
    iter.iter[state], state+1
end

# Dispatch to KnotIterator for nextknotidx and priorknotidx
nextknotidx(iter::KnotRange, x) = nextknotidx(iter.iter, x)
priorknotidx(iter::KnotRange, x) = priorknotidx(iter.iter, x)
