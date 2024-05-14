# Similar to ExtrapDimSpec but for only a single dimension
const ExtrapSpec = Union{BoundaryCondition,Tuple{BoundaryCondition,BoundaryCondition}}

# Union over all BoundaryCondition that yield an infinite sequence of knots
# Note: Extrapolation must result in extrapolated knots, simply evaluating
# outside of the interpolation's bounds is not sufficient
const RepeatKnots = Union{Periodic,Reflect}

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

# Indexing
`KnotIterator` provides limited support for accessing knots via indexing
- `getindex` is provided for `KnotIterator` but does not support Multidimensional
  interpolations (As wrapped by `ProductIterator`) or non-Int indexes.
- A `BoundsError` will be raised if out of bounds and `checkbounds` has been
  implemented for `KnotIterator`

```jldoctest
julia> etp = linear_interpolation([1.0, 1.2, 2.3, 3.0], rand(4); extrapolation_bc=Periodic());

julia> kiter = knots(etp);

julia> kiter[4]
3.0

julia> kiter[36]
24.3
```
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

IteratorSize(::Type{KnotIterator}) = SizeUnknown()
IteratorSize(::Type{KnotIterator{T}}) where {T} = SizeUnknown()

IteratorSize(::Type{KnotIterator{T,ET}}) where {T,ET} = _knot_iter_size(ET)
_knot_iter_size(::Type{<:BoundaryCondition}) = HasLength()
_knot_iter_size(::Type{<:RepeatKnots}) = IsInfinite()
_knot_iter_size(::Type{Tuple{RevBC, FwdBC}}) where {RevBC,FwdBC} = _knot_iter_size(FwdBC)

length(iter::KnotIterator) = _knot_length(iter, IteratorSize(iter))
_knot_length(iter::KnotIterator, ::HasLength) = iter.nknots
size(iter::KnotIterator) = (length(iter),)

IteratorEltype(::Type{KnotIterator}) = EltypeUnknown()
IteratorEltype(::Type{<:KnotIterator{T}}) where {T} = HasEltype()
eltype(::Type{<:KnotIterator{T}}) where {T} = T

"""
    knots(itp::AbstractInterpolation)
    knots(etp::AbstractExtrapolation)

Returns an iterator over knot locations for an AbstractInterpolation or
AbstractExtrapolation.

Iterator will yield scalar values for interpolations over a single dimension,
and tuples of coordinates for higher dimension interpolations. Iteration over
higher dimensions is taken as the product of knots along each dimension.

i.e., Iterator.product(knots on 1st dim, knots on 2nd dim,...)

Extrapolations with Periodic or Reflect boundary conditions, will produce an
infinite sequence of knots.

# Example
```jldoctest
julia> etp = linear_interpolation([1.0, 1.2, 2.3, 3.0], rand(4); extrapolation_bc=Periodic());

julia> Iterators.take(knots(etp), 5) |> collect
5-element Vector{Float64}:
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

# By default start iteration from the first knot, to start at a different index
# use knotsbetween
firstindex(::KnotIterator) = 1
iterate(iter::KnotIterator) = iterate(iter, firstindex(iter))

# Iterate over knots while checkbounds indicates the idx is still valid
function iterate(iter::KnotIterator{T,ET}, idx) where {T, ET}
    if checkbounds(Bool, iter, idx)
        iter[idx], idx+1
    else
        nothing
    end
end

# Split the ExtrapDimSpec into a tuple of left and right boundary conditions
splitExtrapDimSpec(::Type{ET}) where {ET <: BoundaryCondition} = (ET, ET)
splitExtrapDimSpec(::Type{<:Tuple{L,R}}) where {L,R} = (L,R)

# Checks if idx is a valid for the KnotIterator given it's boundary conditions
function checkbounds(::Type{Bool}, iter::KnotIterator{T,ET}, idx::Int) where {T,ET<:ExtrapDimSpec}
    # Get Left/Right Extrapolation Boundary Conditions for the KnotIterator
    left, right = splitExtrapDimSpec(ET)

    # Check Left/Right Boundary Limits -> If Left/Right is a repeating knot, the
    # checks is true as repeated knots generates an infinite sequence of knots
    leftcheck = left <: RepeatKnots ? true : 1 <= idx
    rightcheck = right <: RepeatKnots ? true : idx <= iter.nknots
    leftcheck && rightcheck
end

# Returns the knot given by idx, or raise a BoundsError
# Dispatches to getknotindex(ET, iter, idx) if idx is an extrapolated knot,
# where ET is the boundary conditions for the relevant side
function getindex(iter::KnotIterator{T, ET}, idx::Int) where {T, ET <: ExtrapDimSpec}
    # Get Left/Right Extrapolation Boundary Conditions for the KnotIterator
    left, right = splitExtrapDimSpec(ET)

    # Dispatch to getknotindex if idx is an extrapolated knot, otherwise pull
    # from iter.knots
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

"""
    nextknotidx(iter::KnotIterator, x)

Returns the index of the first knot such that `x < k` or `nothing` if no such
knot exists.

New boundary conditions should define:
```julia
nextknotidx(::Type{<:NewBoundaryCondition}, knots::Vector, x)
```
Where `knots` is `iter.knots` and `NewBoundaryCondition` is the new boundary
conditions. This method is expected to handle values of `x` that are both inbounds
or extrapolated.
"""
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

"""
    priorknotidx(iter::KnotIterator, x)

Returns the index of the last knot such that `k < x` or `nothing` ig not such
knot exists.

New boundary conditions should define
```julia
priorknotidx(::Type{<:NewBoundaryCondition}, knots::Vector, x)
```
Where knots is `iter.knots` and `NewBoundaryCondition` is the new boundary
condition. This method is expected to handle values of `x` that are both inbounds
or extrapolated.
"""
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


# priorknotidx and nextknotidx for non-repeating knots
function nextknotidx(::Type{<:BoundaryCondition}, knots::Vector, x)
    knots[end] <= x ? nothing : findfirst(x .< knots)
end
function priorknotidx(::Type{<:BoundaryCondition}, knots::Vector, x)
    x <= knots[1] ? nothing : findlast(knots .< x)
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

"""
    KnotRange(iter::KnotIterator{T}, start, stop)

Defines an iterator over a range of knots such that `start < k < stop`.

# Fields
- `iter::KnotIterator{T}` Underlying `KnotIterator` providing the knots iterated
- `range::R` Iterator defining the range of knot indices iterated. Where
  `R <: Union{Iterators.Count, UnitRange}`

# Iterator Interface
The following methods defining the Julia's iterator interface have been defined

`Base.IteratorSize` -> Will return one of the following:
- `Base.HasLength` if `range` is of finite length
- `Base.IsInfinite` if `range` is of infinite length
- `Base.SizeUnknown` if the type of `range` is unspecified

`Base.IteratorEltype` -> Returns `Base.EltypeUnknown` if type parameter not
provided, otherwise `Base.HasEltype`

`length` and `size` -> Returns the number of knots to be iterated if
`IteratorSize !== IsInfinite`, otherwise will raise `MethodError`

# Multidimensional Interpolants
Iteration over the knots of a multi-dimensional interpolant is done by wrapping
multiple `KnotRange` iterators within `Iterators.product`.
"""
struct KnotRange{T, R <: Union{Iterators.Count, UnitRange}}
    iter::KnotIterator{T}
    range::R
    KnotRange(iter::KnotIterator{T}, range::R) where {T,R} = new{T,R}(iter, range)
end

function KnotRange(iter::KnotIterator, start, stop)
    # Get the index of the first knot larger than start or firstindex(iter) iff
    # start is missing
    sdx = start === nothing ? firstindex(iter) : nextknotidx(iter, start)

    # Generate knot index range from provided stop / start
    if sdx === nothing
        # Empty Iterator
        range = 1:0
    elseif stop === nothing
        if IteratorSize(iter) === IsInfinite()
            # Stop is missing and iter generates infinite knots
            range = Iterators.countfrom(sdx)
        else
            # Stop is missing and iter generates a finite number of knots
            range = sdx:iter.nknots
        end
    else
        # Stop is provided -> Iterate from start to stop
        edx = priorknotidx(iter, stop)
        range = edx === nothing ? (1:0) : sdx:edx
    end
    KnotRange(iter, range)
end

# Iterator Size is Unknown until we knot range's type
IteratorSize(::Type{KnotRange}) = SizeUnknown()
IteratorSize(::Type{KnotRange{T}}) where {T} = SizeUnknown()
IteratorSize(::Type{KnotRange{T,R}}) where {T, R <: UnitRange} = HasLength()
IteratorSize(::Type{KnotRange{T,R}}) where {T, R <: Iterators.Count} = IsInfinite()

# If type parameter is not provided -> Eltype is Unknown otherwise known
IteratorEltype(::Type{<:KnotRange}) = EltypeUnknown()
IteratorEltype(::Type{<:KnotRange{T}}) where {T} = HasEltype()
eltype(::Type{<:KnotRange{T}}) where {T} = T

# Dispatch length and size to range
length(iter::KnotRange) = length(iter.range)
size(iter::KnotRange) = size(iter.range)

"""
    knotsbetween(iter; start, stop)
    knotsbetween(iter, start, stop)

Iterate over all knots of `iter` such that `start < k < stop`.

`iter` can be an `AbstractInterpolation`, or the output of `knots`
(ie. a `KnotIterator` or `ProductIterator` wrapping `KnotIterator`)

If `start` is not provided, iteration will start from the first knot. An
`ArgumentError` will be raised if both `start` and `stop` are not provided.

If no such knots exists will return a KnotIterator with length 0

# Example
```jldoctest
julia> etp = linear_interpolation([1.0, 1.2, 2.3, 3.0], rand(4); extrapolation_bc=Periodic());

julia> knotsbetween(etp; start=38, stop=42) |> collect
6-element $(Array{typeof(1.0),1}):
 38.3
 39.0
 39.2
 40.3
 41.0
 41.2
```
"""
knotsbetween(itp; start=nothing, stop=nothing) = knotsbetween(itp, start, stop)

# If AbstractInterpolation -> Get Knots -> Then Convert to Knot Range
knotsbetween(itp::AbstractInterpolation, start, stop) = knotsbetween(knots(itp), start, stop)

knotsbetween(::KnotIterator, ::Nothing, ::Nothing) =
    throw(ArgumentError("At least one of `start` or `stop` must be specified"))

knotsbetween(iter::KnotIterator, start, stop) = KnotRange(iter, start, stop)

# Expand Nothing on Start/Stop if tuples
function knotsbetween(iter::Iterators.ProductIterator, start::Union{Nothing, Tuple}, stop::Union{Nothing, Tuple})
    if start === stop === nothing
        throw(ArgumentError("At least one of `start` or `stop` must be specified"))
    else
        kiter = map_ranges(iter.iterators, start, stop)
        Iterators.product(kiter...)
    end
end

# Map iter, start and stop to N KnotRanges, handling start/stop being nothing
function map_ranges(iter, start::Tuple, ::Nothing)
    map((i, s) -> knotsbetween(i, s, nothing), iter, start)
end
function map_ranges(iter, ::Nothing, stop::Tuple)
    map((i, s) -> knotsbetween(i, nothing, s), iter, stop)
end
function map_ranges(iter, start::Tuple, stop::Tuple)
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
