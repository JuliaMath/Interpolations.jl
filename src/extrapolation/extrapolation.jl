struct Extrapolation{T,N,ITPT,IT,ET} <: AbstractExtrapolation{T,N,ITPT,IT}
    itp::ITPT
    et::ET
end

Base.parent(A::Extrapolation) = A.itp
itpflag(etp::Extrapolation) = itpflag(etp.itp)

# DimSpec{Flag} is not enough for extrapolation dispatch, since we allow nested tuples
# However, no tuples should be nested deeper than this; the first level is for different
# schemes in different dimensions, and the second level is for different schemes in
# different directions.
const ExtrapDimSpec = Union{BoundaryCondition,Tuple{Vararg{Union{BoundaryCondition,NTuple{2,BoundaryCondition}}}}}

etptype(::Extrapolation{T,N,ITPT,IT,ET}) where {T,N,ITPT,IT,ET} = ET
etpflag(etp::Extrapolation{T,N,ITPT,IT,ET}) where {T,N,ITPT,IT,ET} = etp.et

"""
`extrapolate(itp, scheme)` adds extrapolation behavior to an interpolation object, according to the provided scheme.

The scheme can take any of these values:

* `Throw` - throws a BoundsError for out-of-bounds indices
* `Flat` - for constant extrapolation, taking the closest in-bounds value
* `Line - linear extrapolation (the wrapped interpolation object must support gradient)
* `Reflect` - reflecting extrapolation (indices must support `mod`)
* `Periodic` - periodic extrapolation (indices must support `mod`)

You can also combine schemes in tuples. For example, the scheme `(Line), Flat())` will use linear extrapolation in the first dimension, and constant in the second.

Finally, you can specify different extrapolation behavior in different direction. `((Line),Flat()), Flat())` will extrapolate linearly in the first dimension if the index is too small, but use constant etrapolation if it is too large, and always use constant extrapolation in the second dimension.
"""
extrapolate(itp::AbstractInterpolation{T,N,IT}, et::ET) where {T,N,IT,ET<:ExtrapDimSpec} =
    Extrapolation{T,N,typeof(itp),IT,ET}(itp, et)

count_interp_dims(::Type{<:Extrapolation{T,N,ITPT}}, n) where {T,N,ITPT} = count_interp_dims(ITPT, n)

@inline function (etp::Extrapolation{T,N})(x::Vararg{Number,N}) where {T,N}
    itp = parent(etp)
    eflag = etpflag(etp)
    xs = inbounds_position(eflag, bounds(itp), x, etp, x)
    extrapolate_value(eflag, skip_flagged_nointerp(itp, x), skip_flagged_nointerp(itp, xs), Tuple(gradient(itp, xs...)), itp(xs...))
end
@inline function (etp::Extrapolation{T,N})(x::Vararg{Union{Number,AbstractVector},N}) where {T,N}
    itp = parent(etp)
    Tret = typeof(lispyprod(zero(T), x...))
    ret = zeros(Tret, shape(x...))
    for (i, y) in zip(eachindex(ret), Iterators.product(x...))
        ret[i] = etp(y...)
    end
    return ret
end

@inline function gradient(etp::AbstractExtrapolation{T,N}, x::Vararg{Number,N}) where {T,N}
    itp = parent(etp)
    if checkbounds(Bool, itp, x...)
        gradient(itp, x...)
    else
        eflag = tcollect(etpflag, etp)
        xs = inbounds_position(eflag, bounds(itp), x, etp, x)
        g = gradient(itp, xs...)
        skipni = t->skip_flagged_nointerp(itp, t)
        SVector(extrapolate_gradient.(skipni(eflag), skipni(x), skipni(xs), Tuple(g)))
    end
end

checkbounds(::Bool, ::AbstractExtrapolation, I...) = true

# The last two arguments are just for error-reporting
function inbounds_position(eflag, bounds, x, etp, xN)
    item = inbounds_index(getfirst(eflag), bounds[1], x[1], etp, xN)
    (item, inbounds_position(getrest(eflag), Base.tail(bounds), Base.tail(x), etp, xN)...)
end
inbounds_position(::Any, ::Tuple{}, ::Tuple{}, etp, xN) = ()

# By default, convert all calls to 2-sided calls
inbounds_index(flag::Flag, bounds, x, etp, xN) = inbounds_index((flag, flag), bounds, x, etp, xN)
# But some one-sided calls can be handled more efficiently that way
function inbounds_index(::Throw, (l,u), x, etp, xN)
    @boundscheck(l <= x <= u || Base.throw_boundserror(etp, xN))
    x
end
inbounds_index(::Periodic, (l,u), x, etp, xN) = periodic(x, l, u)
inbounds_index(::Reflect, (l,u), x, etp, xN) = reflect(x, l, u)

# Left-then-right implementations
function inbounds_index((flagl,flagu)::Tuple{Throw,Flag}, (l,u), x, etp, xN)
    @boundscheck(l <= x || Base.throw_boundserror(etp, xN))
    inbounds_index((nothing,flagu), (l,u), x, etp, xN)
end
function inbounds_index((flagl,flagu)::Tuple{Nothing,Throw}, (l,u), x, etp, xN)
    @boundscheck(x <= u || Base.throw_boundserror(etp, xN))
    x
end

function inbounds_index((flagl,flagu)::Tuple{Union{Flat,Line},Flag}, (l,u), x, etp, xN)
    inbounds_index((nothing,flagu), (l,u), maxp(x,l), etp, xN)
end
function inbounds_index((flagl,flagu)::Tuple{Nothing,Union{Flat,Line}}, (l,u), x, etp, xN)
    minp(x,u)
end

function inbounds_index((flagl,flagu)::Tuple{Periodic,Flag}, (l,u), x, etp, xN)
    inbounds_index((nothing,flagu), (l,u), periodic(x, l, u), etp, xN)
end
function inbounds_index((flagl,flagu)::Tuple{Nothing,Periodic}, (l,u), x, etp, xN)
    periodic(x, l, u)
end

function inbounds_index((flagl,flagu)::Tuple{Reflect,Flag}, (l,u), x, etp, xN)
    inbounds_index((nothing,flagu), (l,u), reflect(x, l, u), etp, xN)
end
function inbounds_index((flagl,flagu)::Tuple{Nothing,Reflect}, (l,u), x, etp, xN)
    reflect(x, l, u)
end

minp(a::T, b::T) where T = min(a, b)
minp(a, b) = min(promote(a, b)...)
maxp(a::T, b::T) where T = max(a, b)
maxp(a, b) = max(promote(a, b)...)

function reflect(y, l, u)
    yr = mod(y - l, 2(u-l)) + l
    return ifelse(yr > u, 2u-yr, yr)
end

periodic(y, l, u) = mod(y-l, u-l) + l


function extrapolate_value(eflag, x, xs, g, val)
    val = extrapolate_axis(getfirst(eflag), x[1], xs[1], g[1], val)
    extrapolate_value(getrest(eflag), Base.tail(x), Base.tail(xs), Base.tail(g), val)
end
extrapolate_value(::Any, ::Tuple{}, ::Tuple{}, ::Tuple{}, val) = val

extrapolate_axis(::Flag, x, xs, g, val) = val
extrapolate_axis(::Line, x, xs, g, val) = val + (x-xs)*g

extrapolate_axis((flagl,flagu)::Tuple{Flag,Flag}, x, xs, g, val) =
    extrapolate_axis((nothing,flagu), x, xs, g, val)
extrapolate_axis((flagl,flagu)::Tuple{Nothing,Flag}, x, xs, g, val) = val

extrapolate_axis((flagl,flagu)::Tuple{Line, Flag}, x, xs, g, val) =
    extrapolate_axis((nothing,flagu), x, xs, g, ifelse(x < xs, val + (x-xs)*g, val))
extrapolate_axis((flagl,flagu)::Tuple{Nothing, Line}, x, xs, g, val) =
    ifelse(x > xs, val + (x-xs)*g, val)

extrapolate_gradient(::Flat, x, xs, g) = ifelse(x==xs, g, zero(g))
extrapolate_gradient(::Flag, x, xs, g) = g

include("filled.jl")
