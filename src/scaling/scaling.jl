export ScaledInterpolation, eachvalue

import Base: iterate

struct ScaledInterpolation{T,N,ITPT,IT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT}
    itp::ITPT
    ranges::RT
end

Base.parent(A::ScaledInterpolation) = A.itp
count_interp_dims(::Type{<:ScaledInterpolation{T,N,ITPT}}, n) where {T,N,ITPT} = count_interp_dims(ITPT, n)

"""
`scale(itp, xs, ys, ...)` scales an existing interpolation object to allow for indexing using other coordinate axes than unit ranges, by wrapping the interpolation object and transforming the indices from the provided axes onto unit ranges upon indexing.

The parameters `xs` etc must be either ranges or linspaces, and there must be one coordinate range/linspace for each dimension of the interpolation object.

For every `NoInterp` dimension of the interpolation object, the range must be exactly `1:size(itp, d)`.
"""
function scale(itp::AbstractInterpolation{T,N,IT}, ranges::Vararg{AbstractRange,N}) where {T,N,IT}
    check_ranges(itpflag(itp), axes(itp), ranges)
    ScaledInterpolation{T,N,typeof(itp),IT,typeof(ranges)}(itp, ranges)
end

function check_ranges(flags, axs, ranges)
    check_range(getfirst(flags), axs[1], ranges[1])
    check_ranges(getrest(flags), Base.tail(axs), Base.tail(ranges))
end
check_ranges(::Any, ::Tuple{}, ::Tuple{}) = nothing

check_range(::NoInterp, ax, r) = ax == r || throw(ArgumentError("The range $r did not equal the corresponding axis of the interpolation object $ax"))
check_range(::Any, ax, r) = length(ax) == length(r) || throw(ArgumentError("The range $r is incommensurate with the corresponding axis $ax"))

# With regards to size and [], ScaledInterpolation behaves like the underlying interpolation object
size(sitp::ScaledInterpolation) = size(sitp.itp)
axes(sitp::ScaledInterpolation) = axes(sitp.itp)

itpflag(sitp::ScaledInterpolation) = itpflag(sitp.itp)

@propagate_inbounds function Base.getindex(sitp::ScaledInterpolation{T,N}, i::Vararg{Int,N}) where {T,N}
    sitp.itp[i...]
end

lbounds(sitp::ScaledInterpolation) = _lbounds(sitp.ranges, itpflag(sitp.itp))
ubounds(sitp::ScaledInterpolation) = _ubounds(sitp.ranges, itpflag(sitp.itp))

boundstep(r::StepRange) = r.step / 2
boundstep(r::UnitRange) = 1//2
"""
Returns *half* the width of one step of the range.

This function is used to calculate the upper and lower bounds of `OnCell` interpolation objects.
""" boundstep

lbound(ax::AbstractRange, ::DegreeBC, ::OnCell) = first(ax) - boundstep(ax)
ubound(ax::AbstractRange, ::DegreeBC, ::OnCell) = last(ax) + boundstep(ax)
lbound(ax::AbstractRange, ::DegreeBC, ::OnGrid) = first(ax)
ubound(ax::AbstractRange, ::DegreeBC, ::OnGrid) = last(ax)

# For (), we scale the evaluation point
function (sitp::ScaledInterpolation{T,N})(xs::Vararg{Number,N}) where {T,N}
    xl = coordslookup(itpflag(sitp.itp), sitp.ranges, xs)
    sitp.itp(xl...)
end
@inline function (sitp::ScaledInterpolation)(x::Vararg{UnexpandedIndexTypes})
    xis = to_indices(sitp, x)
    xis == x && error("evaluation not supported for ScaledInterpolation at positions $x")
    sitp(xis...)
end

(sitp::ScaledInterpolation{T,1}, x::Number, y::Int) where {T} = y == 1 ? sitp(x) : Base.throw_boundserror(sitp, (x, y))

@inline function (itp::ScaledInterpolation{T,N})(x::Vararg{Union{Number,AbstractVector},N}) where {T,N}
    # @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    [itp(i...) for i in Iterators.product(x...)]
end

@inline function coordslookup(flags, ranges, xs)
    item = coordlookup(getfirst(flags), ranges[1], xs[1])
    (item, coordslookup(getrest(flags), Base.tail(ranges), Base.tail(xs))...)
end
coordslookup(::Any, ::Tuple{}, ::Tuple{}) = ()

coordlookup(::NoInterp, r, i) = i
coordlookup(::Flag, r, x) = coordlookup(r, x)

coordlookup(r::UnitRange, x) = x - r.start + oneunit(eltype(r))
# coordlookup(i::Bool, r::AbstractRange, x) = i ? coordlookup(r, x) : convert(typeof(coordlookup(r,x)), x)
coordlookup(r::StepRange, x) = (x - r.start) / r.step + oneunit(eltype(r))

coordlookup(r::StepRangeLen, x) = (x - first(r)) / step(r) + oneunit(eltype(r))
boundstep(r::StepRangeLen) = 0.5*step(r)
rescale_gradient(r::StepRangeLen, g) = g / step(r)

basetype(::Type{ScaledInterpolation{T,N,ITPT,IT,RT}}) where {T,N,ITPT,IT,RT} = ITPT
basetype(sitp::ScaledInterpolation) = basetype(typeof(sitp))


function gradient(sitp::ScaledInterpolation{T,N}, xs::Vararg{Number,N}) where {T,N}
    xl = coordslookup(itpflag(sitp.itp), sitp.ranges, xs)
    g = gradient(sitp.itp, xl...)
    SVector(rescale_gradient_components(itpflag(sitp.itp), sitp.ranges, Tuple(g)))
end

function rescale_gradient_components(flags, ranges, g)
    if getfirst(flags) isa NoInterp
        return rescale_gradient_components(getrest(flags), Base.tail(ranges), g)    # don't consume a coordinate of g
    else
        item = rescale_gradient(ranges[1], g[1])
        return (item, rescale_gradient_components(getrest(flags), Base.tail(ranges), Base.tail(g))...)
    end
end
rescale_gradient_components(flags, ::Tuple{}, ::Tuple{}) = ()

rescale_gradient(r::StepRange, g) = g / r.step
rescale_gradient(r::UnitRange, g) = g

"""
`rescale_gradient(r::AbstractRange)`

Implements the chain rule dy/dx = dy/du * du/dx for use when calculating gradients with scaled interpolation objects.
""" rescale_gradient

### Iteration

struct ScaledIterator{SITPT,CI,WIS}
    sitp::SITPT   # ScaledInterpolation object
    ci::CI        # the CartesianIndices object
    wis::WIS      # WeightedIndex vectors
    breaks1::Vector{Int}   # breaks along dimension 1 where new evaluations must occur
end

Base.IteratorSize(::Type{ScaledIterator{SITPT,CI,WIS}}) where {SITPT,CI<:CartesianIndices{N},WIS} where N = Base.HasShape{N}()
Base.axes(iter::ScaledIterator) = axes(iter.ci)
Base.size(iter::ScaledIterator) = size(iter.ci)

struct ScaledIterState{N,V}
    cistate::CartesianIndex{N}
    ibreak::Int
    cached_evaluations::NTuple{N,V}
end

function eachvalue(sitp::ScaledInterpolation{T,N}) where {T,N}
    itps = tcollect(itpflag, sitp.itp)
    newaxes = map(r->Base.Slice(ceil(Int, first(r)):floor(Int, last(r))), sitp.ranges)
    wis = dimension_wis(value_weights, itps, axes(sitp.itp), newaxes, sitp.ranges)
    wis1 = wis[1]
    i1 = first(axes(wis1, 1))
    breaks1 = [i1]
    for i in Iterators.drop(axes(wis1, 1), 1)
        if indexes(wis1[i]) != indexes(wis1[i-1])
            push!(breaks1, i)
        end
    end
    push!(breaks1, last(axes(wis1, 1))+1)
    ScaledIterator(sitp, CartesianIndices(newaxes), wis, breaks1)
end

function dimension_wis(f::F, itps, axs, newaxes, ranges) where F
    itpflag, ax, nax, r = itps[1], axs[1], newaxes[1], ranges[1]
    function makewi(x)
        pos, coefs = weightedindex_parts((f,), itpflag, ax, coordlookup(r, x))
        maybe_weightedindex(pos, coefs[1])
    end
    (makewi.(nax), dimension_wis(f, Base.tail(itps), Base.tail(axs), Base.tail(newaxes), Base.tail(ranges))...)
end
dimension_wis(f, ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Tuple{}) = ()

function Base.iterate(iter::ScaledIterator)
    ret = iterate(iter.ci)
    ret === nothing && return nothing
    item, cistate = ret
    wis = getindex.(iter.wis, Tuple(item))
    ces = cache_evaluations(iter.sitp.itp.coefs, indexes(wis[1]), weights(wis[1]), Base.tail(wis))
    return _reduce(+, weights(wis[1]).*ces), ScaledIterState(cistate, first(iter.breaks1), ces)
end

function Base.iterate(iter::ScaledIterator, state)
    ret = iterate(iter.ci, state.cistate)
    ret === nothing && return nothing
    item, cistate = ret
    i1 = item[1]
    isnext1 = i1 == state.cistate[1]+1
    if isnext1 && i1 < iter.breaks1[state.ibreak+1]
        # We can use the previously cached values
        wis1 = iter.wis[1][i1]
        return _reduce(+, weights(wis1).*state.cached_evaluations), ScaledIterState(cistate, state.ibreak, state.cached_evaluations)
    end
    # Re-evaluate. We're being a bit lazy here: in some cases, some of the cached values could be reused
    wis = getindex.(iter.wis, Tuple(item))
    ces = cache_evaluations(iter.sitp.itp.coefs, indexes(wis[1]), weights(wis[1]), Base.tail(wis))
    return _reduce(+, weights(wis[1]).*ces), ScaledIterState(cistate, isnext1 ? state.ibreak+1 : first(iter.breaks1), ces)
end

_reduce(op, list) = op(list[1], _reduce(op, Base.tail(list)))
_reduce(op, list::Tuple{Number}) = list[1]
_reduce(op, list::Tuple{}) = error("cannot reduce an empty list")

# We use weights only as a ruler to determine when we are done
cache_evaluations(coefs, i::Int, weights, rest) = (coefs[i, rest...], cache_evaluations(coefs, i+1, Base.tail(weights), rest)...)
cache_evaluations(coefs, indexes, weights, rest) = (coefs[indexes[1], rest...], cache_evaluations(coefs, Base.tail(indexes), Base.tail(weights), rest)...)
cache_evaluations(coefs, ::Int, ::Tuple{}, rest) = ()
cache_evaluations(coefs, ::Any, ::Tuple{}, rest) = ()

ssize(sitp::ScaledInterpolation{T,N}) where {T,N} = map(r->round(Int, last(r)-first(r)+1), sitp.ranges)::NTuple{N,Int}
