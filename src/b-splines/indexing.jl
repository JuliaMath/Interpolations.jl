### Primary evaluation (indexing) entry points

@inline function (itp::BSplineInterpolation{T,N})(x::Vararg{Number,N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    expand_value(itp, x)
end
@propagate_inbounds function (itp::BSplineInterpolation{T,N})(x::Vararg{Number,M}) where {T,M,N}
    inds, trailing = split_trailing(itp, x)
    @boundscheck (check1(trailing) || Base.throw_boundserror(itp, x))
    @assert length(inds) == N
    itp(inds...)
end

@inline function gradient(itp::BSplineInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    @boundscheck checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x)
    expand_gradient(itp, x)
end
@inline function gradient!(dest, itp::BSplineInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    @boundscheck checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x)
    expand_gradient!(dest, itp, x)
end

@inline function hessian(itp::BSplineInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    @boundscheck checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x)
    expand_hessian(itp, x)
end
@inline function hessian!(dest, itp::BSplineInterpolation{T,N}, x::Vararg{Number,N}) where {T,N}
    @boundscheck checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x)
    expand_hessian(itp, x)
end

checkbounds(::Type{Bool}, itp::AbstractInterpolation, x::Vararg{Number,N}) where N =
    checklubounds(lbounds(itp), ubounds(itp), x)

# Leftovers from AbstractInterpolation
@inline function (itp::BSplineInterpolation)(x::Vararg{UnexpandedIndexTypes})
    itp(to_indices(itp, x)...)
end
@inline function (itp::BSplineInterpolation)(x::Vararg{ExpandedIndexTypes})
    itp.(Iterators.product(x...))
end

"""
    val = expand_value(itp, x)

Interpolate `itp` at `x`.
"""
function expand_value(itp::AbstractInterpolation, x::Tuple)
    coefs = coefficients(itp)
    degree = interpdegree(itp)
    ixs, rxs = splitgrouped(expand_indices_resid(degree, axes(itp), x))
    cxs = expand_weights(value_weights, degree, rxs)
    expand(coefs, cxs, ixs)
end

"""
    g = expand_gradient(itp, x)

Calculate the interpolated gradient of `itp` at `x`.
"""
function expand_gradient(itp::AbstractInterpolation, x::Tuple)
    coefs = coefficients(itp)
    degree = interpdegree(itp)
    ixs, rxs = splitgrouped(expand_indices_resid(degree, axes(itp), x))
    cxs = expand_weights(value_weights, degree, rxs)
    gxs = expand_weights(gradient_weights, degree, rxs)
    expand(coefs, (cxs, gxs), ixs)
end

function expand_gradient!(dest, itp::AbstractInterpolation, x::Tuple)
    coefs = coefficients(itp)
    degree = interpdegree(itp)
    ixs, rxs = splitgrouped(expand_indices_resid(degree, axes(itp), x))
    cxs = expand_weights(value_weights, degree, rxs)
    gxs = expand_weights(gradient_weights, degree, rxs)
    expand!(dest, coefs, (cxs, gxs), ixs)
end

"""
    H = expand_hessian(itp, x)

Calculate the interpolated hessian of `itp` at `x`.
"""
function expand_hessian(itp::AbstractInterpolation, x::Tuple)
    coefs = coefficients(itp)
    degree = interpdegree(itp)
    ixs, rxs = splitgrouped(expand_indices_resid(degree, axes(itp), x))
    cxs = expand_weights(value_weights, degree, rxs)
    gxs = expand_weights(gradient_weights, degree, rxs)
    hxs = expand_weights(hessian_weights, degree, rxs)
    expand(coefs, (cxs, gxs, hxs), ixs)
end

"""
    Weights{N}

A type alias for an `N`-dimensional tuple of weights or indexes for interpolation.

# Example

If you are performing linear interpolation (degree 1) in three dimensions at the point
`x = [5.1, 8.2, 4.3]`, then the floating-point residuals are `[0.1, 0.2, 0.3]`.
For value interpoations, the corresponding `Weights` would be

    ((0.9,0.1),  (0.8,0.2),  (0.7,0.3))    # (dim1,  dim2,  dim3)

Note each "inner" tuple, for value interpolation, sums to 1. For gradient Weights,
each inner tuple would be `(-1, 1)`, and for hessian Weights each would be `(0, 0)`.

The same structure can be used for the integer indexes at which each coefficient is evaluated.
For the example above (with `x = [5.1, 8.2, 4.3]`), the indexes would be

    ((5,6), (8,9), (4,5))

corresponding to the integer pairs that bracket each coordinate.

When performing mixed interpolation (e.g., `Linear` along dimension 1 and `Cubic` along dimension 2),
the inner tuples may not all be of the same length.
"""
const Weights{N} = NTuple{N,Tuple{Vararg{<:Number}}}

"""
    Indexes{N}

The same as [`Weights`](@ref) for the integer-values used as evaluation locations for the
coefficients array.
"""
const Indexes{N} = NTuple{N,Tuple{Vararg{<:Integer}}}

"""
    val = expand(coefs, vweights::Weights, ixs::Indexes)
    g   = expand(coefs, (vweights, gweights), ixs)
    H   = expand(coefs, (vweights, gweights, hweights), ixs)

Calculate the value, gradient, or hessian of a separable AbstractInterpolation object.
This function works recursively, processing one index at a time and calling itself as

    ret = expand(coefs, <weights>, ixs, iexpanded...)

(The `iexpanded` form is not intended to be called directly by users.) For example,
for two-dimensional linear interpolation at a point `x = [5.1, 8.2]` the corresponding
top-level call would be

    #              weights                  ixs
    expand(coefs,  ((0.9,0.1), (0.8,0.2)),  ((5,6), (8,9))

After one round of recursion this becomes

    #                  weights        ixs        iexpanded
    0.9*expand(coefs,  ((0.8,0.2),),  ((8,9),),  5) +
    0.1*expand(coefs,  ((0.8,0.2),),  ((8,9),),  6)

(The first dimension has been processed.) After another round, this becomes

    #                      wts ixs iexpanded
    0.9*(0.8*expand(coefs, (), (), 5, 8) + 0.2*expand(coefs, (), (), 5, 9)) +
    0.1*(0.8*expand(coefs, (), (), 6, 8) + 0.2*expand(coefs, (), (), 6, 9))

Now that the weights and `ixs` are empty and all indices are in `iexpanded`,
it finally resolves to

    0.9*(0.8*coefs[5, 8] + 0.2*coefs[5, 9]) +
    0.1*(0.8*coefs[6, 8] + 0.2*coefs[6, 9])

which is the expression for bilinear interpolation at the given `x`.

For calculating the components of the gradient and hessian, individual dimensions of
`gweights` and/or `hweights` will be substituted into the appropriate slot. For example,
in three dimensions

    g[1] = expand(coefs, (gweights[1], vweights[2], vweights[3]), ixs)
    g[2] = expand(coefs, (vweights[1], gweights[2], vweights[3]), ixs)
    g[3] = expand(coefs, (vweights[1], vweights[2], gweights[3]), ixs)
"""
function expand(coefs::AbstractArray, vweights::Weights, ixs::Indexes, iexpanded::Vararg{Integer,M}) where {M}
    w1, wrest = vweights[1], Base.tail(vweights)
    ix1, ixrest = ixs[1], Base.tail(ixs)
    _expand1(coefs, w1, ix1, wrest, ixrest, iexpanded)
end
function expand(coefs::AbstractArray{T,N}, vweights::Tuple{}, ixs::Tuple{}, iexpanded::Vararg{Integer,N}) where {T,N}
    @inbounds coefs[iexpanded...]  # @inbounds is safe because we checked in the original call
end

const HasNoInterp{N} = NTuple{N,Tuple{Vararg{<:Union{Number,NoInterp}}}}
expand(coefs::AbstractArray, vweights::HasNoInterp, ixs::Indexes, iexpanded::Vararg{Integer,M}) where {M} = NoInterp()

# _expand1 handles the expansion of a single dimension weight list (of length L)
@inline _expand1(coefs, w1, ix1, wrest, ixrest, iexpanded) =
    w1[1] * expand(coefs, wrest, ixrest, iexpanded..., ix1[1]) +
    _expand1(coefs, Base.tail(w1), Base.tail(ix1), wrest, ixrest, iexpanded)
@inline _expand1(coefs, w1::Tuple{Number}, ix1::Tuple{Integer}, wrest, ixrest, iexpanded) =
    w1[1] * expand(coefs, wrest, ixrest, iexpanded..., ix1[1])

# Expansion of the gradient
function expand(coefs, (vweights, gweights)::Tuple{HasNoInterp{N},HasNoInterp{N}}, ixs::Indexes{N}) where N
    # We swap in one gradient dimension per call to expand
    SVector(skip_nointerp(ntuple(d->expand(coefs, substitute(vweights, d, gweights), ixs), Val(N))...))
end
function expand!(dest, coefs, (vweights, gweights)::Tuple{HasNoInterp{N},HasNoInterp{N}}, ixs::Indexes{N}) where N
    # We swap in one gradient dimension per call to expand
    i = 0
    for d = 1:N
        w = substitute(vweights, d, gweights)
        w isa Weights || continue   # if this isn't true, it must have a NoInterp in it
        dest[i+=1] = expand(coefs, w, ixs)
    end
    dest
end

# Expansion of the hessian
# To handle the immutability of SMatrix we build static methods that visit just the entries we need,
# which due to symmetry is just the upper triangular part
ntuple_sym(f, ::Val{0}) = ()
ntuple_sym(f, ::Val{1}) = (f(1,1),)
ntuple_sym(f, ::Val{2}) = (f(1,1), f(1,2), f(2,2))
ntuple_sym(f, ::Val{3}) = (f(1,1), f(1,2), f(2,2), f(1,3), f(2,3), f(3,3))
ntuple_sym(f, ::Val{4}) = (f(1,1), f(1,2), f(2,2), f(1,3), f(2,3), f(3,3), f(1,4), f(2,4), f(3,4), f(4,4))
@inline function ntuple_sym(f, ::Val{N}) where N
    (ntuple_sym(f, Val(N-1))..., ntuple(i->f(i,N), Val(N))...)
end

sym2dense(t::Tuple{})              = t
sym2dense(t::NTuple{1,T}) where T  = t
sym2dense(t::NTuple{3,T}) where T  = (t[1], t[2], t[2], t[3])
sym2dense(t::NTuple{6,T}) where T  = (t[1], t[2], t[4], t[2], t[3], t[5], t[4], t[5], t[6])
sym2dense(t::NTuple{10,T}) where T = (t[1], t[2], t[4], t[7], t[2], t[3], t[5], t[8], t[4], t[5], t[6], t[9], t[7], t[8], t[9], t[10])
function sym2dense(t::NTuple{L,T}) where {L,T}
    # Warning: non-inferrable unless we make this @generated.
    # Above 4 dims one might anyway prefer an Array, and use hessian!
    N = ceil(Int, sqrt(2*L))
    @assert (N*(N+1))÷2 == L
    a = Vector{T}(undef, N*N)
    idx = 0
    for j = 1:N, i=1:N
        iu, ju = ifelse(i>=j, (j, i), (i, j))  # index in the upper triangular
        k = (ju*(ju+1))÷2 + iu
        a[idx+=1] = t[k]
    end
    tuple(a...)
end

squarematrix(t::NTuple{1,T}) where T  = SMatrix{1,1,T}(t)
squarematrix(t::NTuple{4,T}) where T  = SMatrix{2,2,T}(t)
squarematrix(t::NTuple{9,T}) where T  = SMatrix{3,3,T}(t)
squarematrix(t::NTuple{16,T}) where T = SMatrix{4,4,T}(t)
function squarematrix(t::NTuple{L,T}) where {L,T}
    # Warning: non-inferrable unless we make this @generated.
    # Above 4 dims one might anyway prefer an Array, and use hessian!
    N = floor(Int, sqrt(L))
    @assert N*N == L
    SMatrix{N,N,T}(t)
end

cumweights(w) = _cumweights(0, w...)
_cumweights(c, w1, w...) = (c+1, _cumweights(c+1, w...)...)
_cumweights(c, ::NoInterp, w...) = (c, _cumweights(c, w...)...)
_cumweights(c) = ()

function expand(coefs, (vweights, gweights, hweights)::NTuple{3,HasNoInterp{N}}, ixs::Indexes{N}) where N
    coefs = ntuple_sym((i,j)->expand(coefs, substitute(vweights, i, j, gweights, hweights), ixs), Val(N))
    squarematrix(sym2dense(skip_nointerp(coefs...)))
end

function expand!(dest, coefs, (vweights, gweights, hweights)::NTuple{3,HasNoInterp{N}}, ixs::Indexes{N}) where N
    # The Hessian is nominally N × N, but if there are K NoInterp dims then it's N-K × N-K
    indlookup = cumweights(hweights)   # for d in 1:N, indlookup[d] returns the appropriate index in 1:N-K
    for d2 = 1:N, d1 = 1:d2
        w = substitute(vweights, d1, d2, gweights, hweights)
        w isa Weights || continue   # if this isn't true, it must have a NoInterp in it
        i, j = indlookup[d1], indlookup[d2]
        dest[i, j] = dest[j, i] = expand(coefs, w, ixs)
    end
    dest
end

function expand_indices_resid(degree, axs, x)
    item = expand_index_resid(getfirst(degree), axs[1], x[1])
    (item, expand_indices_resid(getrest(degree), Base.tail(axs), Base.tail(x))...)
end
expand_indices_resid(degree, ::Tuple{}, ::Tuple{}) = ()

function expand_index_resid(degree, ax, x::Number)
    ix, δx = base_rem(degree, ax, x)
    expand_index(degree, ix, ax, δx), δx
end

expand_weights(f, degree::Union{Degree,NoInterp}, ixs) =
    (f(degree, ixs[1]), expand_weights(f, degree, Base.tail(ixs))...)
expand_weights(f, degree::Union{Degree,NoInterp}, ::Tuple{}) = ()

expand_weights(f, degree::Tuple{Vararg{Union{Degree,NoInterp},N}}, ixs::NTuple{N,Number}) where N =
    f.(degree, ixs)

# expand_indices(degree::Union{Degree,NoInterp}, ixs, axs, δxs) =
#     (expand_index(degree, ixs[1], axs[1], δxs[1]), expand_indices(degree, Base.tail(ixs), Base.tail(axs), Base.tail(δxs))...)
# expand_indices(degree::Union{Degree,NoInterp}, ::Tuple{}, ::Tuple{}, ::Tuple{}) = ()

# expand_indices(degree::Tuple{Vararg{Union{Degree,NoInterp},N}}, ixs::NTuple{N,Number}, axs::NTuple{N,Tuple{Real,Real}}, δxs::NTuple{N,Number}) where N =
#     expand_index.(degree, ixs, axs, δxs)

# expand_index(degree, ixs, bounds::Tuple{Real,Real}, δxs) = expand_index(degree, ixs, axfrombounds(bounds), δxs)

checklubounds(ls, us, xs) = _checklubounds(true, ls, us, xs)
_checklubounds(tf::Bool, ls, us, xs) = _checklubounds(tf & (ls[1] <= xs[1] <= us[1]),
                                                      Base.tail(ls), Base.tail(us), Base.tail(xs))
_checklubounds(tf::Bool, ::Tuple{}, ::Tuple{}, ::Tuple{}) = tf


# there is a Heisenbug, when Base.promote_op is inlined into getindex_return_type
# thats why we use this @noinline fence
@noinline _promote_mul(a,b) = Base.promote_op(*, a, b)

@noinline function getindex_return_type(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}}, argtypes::Tuple) where {T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}
    reduce(_promote_mul, eltype(TCoefs), argtypes)
end

function getindex_return_type(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}}, ::Type{I}) where {T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad,I}
    _promote_mul(eltype(TCoefs), I)
end

# This handles round-towards-the-middle for points on half-integer edges
roundbounds(x::Integer, bounds) = x
function roundbounds(x, bounds)
    l, u = first(bounds), last(bounds)
    h = half(x)
    xh = x+h
    ifelse(x < u+half(u), floor(xh), ceil(xh)-1)
end

floorbounds(x::Integer, ax) = x
function floorbounds(x, ax)
    l = first(ax)
    h = half(x)
    ifelse(x < l, floor(x+h), floor(x+zero(h)))
end

half(x) = oneunit(x)/2
