# Indexing at a point
@inline function (itp::GriddedInterpolation{T,N})(x::Vararg{Number,N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    wis = weightedindexes((value_weights,), itpinfo(itp)..., x)
    coefficients(itp)[wis...]
end
@propagate_inbounds function (itp::GriddedInterpolation{T,N})(x::Vararg{Number,M}) where {T,M,N}
    inds, trailing = split_trailing(itp, x)
    @boundscheck (check1(trailing) || Base.throw_boundserror(itp, x))
    @assert length(inds) == N
    itp(inds...)
end
@inline function (itp::GriddedInterpolation)(x::Vararg{UnexpandedIndexTypes})
    xis = to_indices(itp, x)
    xis == x && error("evaluation not supported for GriddedInterpolation at positions $x")
    itp(xis...)
end

@inline function gradient(itp::GriddedInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    wis = weightedindexes((value_weights, gradient_weights), itpinfo(itp)..., x)
    SVector(map(inds->coefficients(itp)[inds...], wis))
end

itpinfo(itp::GriddedInterpolation) = (tcollect(itpflag, itp), itp.knots)

roundbounds(x::Integer, knotvec::AbstractVector) = gridded_roundbounds(x, knotvec)
roundbounds(x::Number, knotvec::AbstractVector) = gridded_roundbounds(x, knotvec)
function gridded_roundbounds(x, knotvec::AbstractVector)
    i = find_knot_index(knotvec, x)
    iclamp = max(i, first(axes1(knotvec)))
    inext = min(iclamp+1, last(axes1(knotvec)))
    ifelse(i < iclamp, i+1, ifelse(x - knotvec[iclamp] < knotvec[inext] - x, i, inext))
end

floorbounds(x::Integer, knotvec::AbstractVector) = gridded_floorbounds(x, knotvec)
floorbounds(x, knotvec::AbstractVector) = gridded_floorbounds(x, knotvec)
function gridded_floorbounds(x, knotvec::AbstractVector)
    i = find_knot_index(knotvec, x)
    max(i, first(axes1(knotvec)))
end

@inline find_knot_index(knotv, x) = searchsortedfirst(knotv, x, Base.Order.ForwardOrdering()) - 1
@inline find_knot_index(knotv, x::AbstractVector) = searchsortedfirst_vec(knotv, x) .- 1

@inline function weightedindex_parts(fs::F, mode::Gridded, knotvec::AbstractVector, x) where F
    i = find_knot_index(knotvec, x)
    weightedindex_parts2(fs, mode, knotvec, x, i)
end

@inline function weightedindex_parts2(fs::F, mode::Gridded, knotvec::AbstractVector, x, i) where F
    ax1 = axes1(knotvec)
    iclamp = clamp(i, first(ax1), last(ax1)-1)
    weightedindex(fs, degree(mode), knotvec, x, iclamp)
end

function weightedindex(fs::F, deg::Constant, knotvec, x, iclamp) where F
    pos, δx = positions(deg, knotvec,  x)
    (position=pos, coefs=fmap(fs, deg, δx))
end
function weightedindex(fs::F, deg::Degree, knotvec, x, iclamp) where F
    @inbounds l, u = knotvec[iclamp], knotvec[iclamp+1]
    δx = ratio(x - l, u - l)
    (position=iclamp, coefs=rescale_gridded(fs, fmap(fs, deg, δx), u-l))
end

rescale_gridded(fs::F, coefs, Δx) where F =
    (rescale_gridded(fs[1], coefs[1], Δx), rescale_gridded(Base.tail(fs), Base.tail(coefs), Δx)...)
rescale_gridded(::Tuple{}, ::Tuple{}, Δx) = ()
rescale_gridded(::typeof(value_weights), coefs, Δx) = coefs
rescale_gridded(::typeof(gradient_weights), coefs, Δx) = coefs./Δx
rescale_gridded(::typeof(hessian_weights), coefs, Δx) = coefs./Δx.^2

@inline function (itp::GriddedInterpolation{T,N})(x::Vararg{Union{Number,AbstractVector},N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    itps = tcollect(itpflag, itp)
    if x[1] isa AbstractVector
        wis = dimension_wis_vec(value_weights, itps, itp.knots, x)
    else
        wis = dimension_wis(value_weights, itps, itp.knots, x)
    end
    coefs = coefficients(itp)
    ret = [coefs[i...] for i in Iterators.product(wis...)]
    reshape(ret, shape(wis...))
end

function dimension_wis(f::F, itps, knots, xs) where F
    itpflag, knotvec, x = itps[1], knots[1], xs[1]
    function makewi(y)
        pos, coefs = weightedindex_parts((f,), itpflag, knotvec, y)
        maybe_weightedindex(pos, coefs[1])
    end
    (makewi.(x), dimension_wis(f, Base.tail(itps), Base.tail(knots), Base.tail(xs))...)
end

function dimension_wis_vec(f::F, itps, knots, xs) where F
    itpflag, knotvec, x = itps[1], knots[1], xs[1]
    ivec = find_knot_index(knotvec, x)
    function makewi(y, i)
        pos, coefs = weightedindex_parts2((f,), itpflag, knotvec, y, i)
        maybe_weightedindex(pos, coefs[1])
    end
    (makewi.(x, ivec), dimension_wis(f, Base.tail(itps), Base.tail(knots), Base.tail(xs))...)
end

function dimension_wis(f::F, itps::Tuple{NoInterp,Vararg{Any}}, knots, xs) where F
    (Int.(xs[1]), dimension_wis(f, Base.tail(itps), Base.tail(knots), Base.tail(xs))...)
end
dimension_wis(f, ::Tuple{}, ::Tuple{}, ::Tuple{}) = ()


function getindex_return_type(::Type{GriddedInterpolation{T,N,TCoefs,IT,K}}, argtypes) where {T,N,TCoefs,IT<:DimSpec{Gridded},K}
    Tret = TCoefs
    for a in argtypes
        Tret = Base.promote_op(*, Tret, a)
    end
    Tret
end

Base.@propagate_inbounds function searchsortedfirst_exp_left(v, xx, lo, hi)
    for i in 0:4
        ind = lo + i
        ind > hi && return ind
        xx <= v[ind] && return ind
    end
    n = 3
    tn2 = 2^n
    tn2m1 = 2^(n-1)
    ind = lo + tn2
    while ind <= hi
        xx <= v[ind] && return searchsortedfirst(v, xx, lo + tn2 - tn2m1, ind, Base.Order.Forward)
        tn2 *= 2
        tn2m1 *= 2
        ind = lo + tn2
    end
    return searchsortedfirst(v, xx, lo + tn2 - tn2m1, hi, Base.Order.Forward)
end

function searchsortedfirst_vec(v::AbstractVector, x::AbstractVector)
    issorted(x) || return searchsortedfirst.(Ref(v), x)
    out = zeros(Int, length(x))
    lo = 1
    hi = length(v)
    @inbounds for i in 1:length(x)
        xx = x[i]
        y = searchsortedfirst_exp_left(v, xx, lo, hi)
        out[i] = y
        lo = min(y, hi)
    end
    return out
end
