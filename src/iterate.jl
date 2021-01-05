
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

function knots(itp::AbstractInterpolation)
    I = Iterators.product(getknots(itp)...)
    ndims(itp) == 1 ? Iterators.flatten(I) : I
end

function knots(etp::AbstractExtrapolation)
    k = getknots(etp)
    bc = etpflag(etp)
    iter = KnotIterator(k, bc)
    length(iter) == 1 ? only(iter) : Iterators.ProductIterator(iter)
end

Base.iterate(iter::KnotIterator{T}) where {T} = iterate(iter, (1, zero(T)))

function Base.iterate(iter::KnotIterator, state)
    curidx, offset = state[1], state[2]
    isnothing(curidx) && return nothing
    idx = inbounds_index(iter.bc, (1, iter.nknots), curidx, iter, state)
    iter.knots[idx] + offset, nextstate(iter, state)
end


function nextstate(iter::KnotIterator, state)
    idx, offset = state
    nextidx = idx < iter.nknots ? idx + 1 : nothing
    (nextidx, offset)
end
function nextstate(iter::KnotIterator{T,ET}, state) where {T,ET<:RepeatKnots}
    idx, offset = state
    if idx + 1 < iter.nknots
        return idx + 1, offset
    else
        knotrange = iter.knots[end] - iter.knots[1]
        return 1, offset + knotrange
    end
end
