@inline sqr(x) = x*x
@inline cub(x) = x*x*x

modrange(x, r::AbstractUnitRange) = mod(x-first(r), length(r)) + first(r)
modrange(x, (l, u)::Tuple{Real,Real}) = mod(x-l, u-l+1) + l

fmap(fs, x...) = _fmap(x, fs...)
@inline _fmap(x, f, fs...) = (f(x...), _fmap(x, fs...)...)
@inline _fmap(x) = ()

split_flag(f::Flag) = f, f
split_flag(t::Tuple) = t[1], Base.tail(t)

getfirst(f::Flag) = f
getfirst(t::Tuple) = t[1]
getrest(f::Flag) = f
getrest(t::Tuple) = Base.tail(t)

tcollect(f, itp::AbstractInterpolation{T,N}) where {T,N} = _tcollect(ntuple(d->true, Val(N)), f(itp))
@inline _tcollect(ruler, prop) = (getfirst(prop), _tcollect(Base.tail(ruler), getrest(prop))...)
_tcollect(::Tuple{}, prop) = ()

split_trailing(::AbstractArray{T,N}, x) where {T,N} = Base.IteratorsMD.split(x, Val(N))
check1(args) = _check1(true, args...)
@inline _check1(tf, a, args...) = _check1(tf & (a == 1), args...)
_check1(tf) = tf

# These are not inferrable for mixed-type tuples, so when that's important use `getfirst`
# and `getrest` instead.
iextract(f::Flag, d) = f
iextract(t::Tuple, d) = t[d]

splitgrouped(prs::Tuple{Vararg{NTuple{2,Any}}}) = first.(prs), last.(prs)
splitgrouped(prs::Tuple{Vararg{NTuple{3,Any}}}) = first.(prs), middle.(prs), last.(prs)
middle(t::Tuple{Any,Any,Any}) = t[2]

fast_trunc(::Type{Int}, x) = unsafe_trunc(Int, x)
fast_trunc(::Type{Int}, x::Rational) = x.num ÷ x.den

# Slot-substitution guided by a `ruler` tuple. Substitution occurs when `default` has the same
# length as `ruler`.
@inline substitute_ruled(default, ruler, subst) = (default[1], substitute_ruled(Base.tail(default), ruler, Base.tail(subst))...)
@inline substitute_ruled(default::NTuple{N,Any}, ruler::NTuple{N,Any}, subst) where N =
    (subst[1], substitute_ruled(Base.tail(default), ruler, Base.tail(subst))...)
substitute_ruled(default::Tuple{}, ruler::NTuple{N,Any}, subst) where N = ()

@inline skip_nointerp(x, rest...) = (x, skip_nointerp(rest...)...)
@inline skip_nointerp(::NoInterp, rest...) = skip_nointerp(rest...)
skip_nointerp() = ()

@inline sumvals(val, δval, args...) = sumvals(val+δval, args...)
@inline sumvals(val, ::Nothing, args...) = sumvals(val, args...)
sumvals(val) = val

@inline promote_typeof(a, b, args...) = _promote_typeof(promote_type(typeof(a), typeof(b)), args...)
@inline _promote_typeof(::Type{T}, a, args...) where T = _promote_typeof(promote_type(T, typeof(a)), args...)
_promote_typeof(::Type{T}) where T = T
