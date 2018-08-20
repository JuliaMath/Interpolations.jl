@inline sqr(x) = x*x
@inline cub(x) = x*x*x

modrange(x, r::AbstractUnitRange) = mod(x-first(r), length(r)) + first(r)
modrange(x, (l, u)::Tuple{Real,Real}) = mod(x-l, u-l+1) + l

split_flag(f::Flag) = f, f
split_flag(t::Tuple) = t[1], Base.tail(t)

getfirst(f::Flag) = f
getfirst(t::Tuple) = t[1]
getrest(f::Flag) = f
getrest(t::Tuple) = Base.tail(t)

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

# Substitution for gradient components
function substitute(default::NTuple{N,Any}, d::Integer, subst::NTuple{N,Any}) where N
    ntuple(i->ifelse(i==d, subst[i], default[i]), Val(N))
end
function substitute(default::NTuple{N,Any}, d::Integer, val) where N
    ntuple(i->ifelse(i==d, val, default[i]), Val(N))
end

# Substitution for hessian components
function substitute(default::NTuple{N,Any}, d1::Integer, d2::Integer, subst1::NTuple{N,Any}, subst2::NTuple{N,Any}) where N
    ntuple(i->ifelse(i==d1==d2, subst2[i], ifelse(i==d1, subst1[i], ifelse(i==d2, subst1[i], default[i]))), Val(N))
end

@inline skip_nointerp(x, rest...) = (x, skip_nointerp(rest...)...)
@inline skip_nointerp(::NoInterp, rest...) = skip_nointerp(rest...)
skip_nointerp() = ()

@inline sumvals(val, δval, args...) = sumvals(val+δval, args...)
@inline sumvals(val, ::Nothing, args...) = sumvals(val, args...)
sumvals(val) = val

@inline promote_typeof(a, b, args...) = _promote_typeof(promote_type(typeof(a), typeof(b)), args...)
@inline _promote_typeof(::Type{T}, a, args...) where T = _promote_typeof(promote_type(T, typeof(a)), args...)
_promote_typeof(::Type{T}) where T = T
