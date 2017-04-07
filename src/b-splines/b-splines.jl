export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic,
    Cubic

@compat abstract type Degree{N} <: Flag end

immutable BSpline{D<:Degree} <: InterpolationType end
BSpline{D<:Degree}(::D) = BSpline{D}()

bsplinetype{D<:Degree}(::Type{BSpline{D}}) = D

immutable BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad} <: AbstractInterpolation{T,N,IT,GT}
    coefs::TCoefs
end
function BSplineInterpolation{N,Tel,TWeights<:Real,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad}(::Type{TWeights}, A::AbstractArray{Tel,N}, ::IT, ::GT, ::Val{pad})
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(typeof(A)) || warn("For performance reasons, consider using an array of a concrete type (typeof(A) == $(typeof(A)))")

    c = one(TWeights)
    for _ in 2:N
        c *= c
    end
    T = typeof(c * one(Tel))

    BSplineInterpolation{T,N,typeof(A),IT,GT,pad}(A)
end

# Utilities for working either with scalars or tuples/tuple-types
iextract{T<:BSpline}(::Type{T}, d) = T
iextract(t, d) = t.parameters[d]
padding{T,N,TCoefs,IT,GT,pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}) = pad
padding(itp::AbstractInterpolation) = padding(typeof(itp))
padextract(pad::Integer, d) = pad
padextract(pad::Tuple{Vararg{Integer}}, d) = pad[d]

lbound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d) =
    first(indices(itp, d))
ubound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d) =
    last(indices(itp, d))
lbound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d) =
    first(indices(itp, d)) - 0.5
ubound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d) =
    last(indices(itp, d))+0.5

count_interp_dims{T,N,TCoefs,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType},pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}, n) = count_interp_dims(IT, n)

function size{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d)
    d <= N ? size(itp.coefs, d) - 2*padextract(pad, d) : 1
end

@inline indices{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}) =
    map_repeat(indices_removepad, indices(itp.coefs), pad)

function indices{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d)
    d <= N ? indices_removepad(indices(itp.coefs, d), padextract(pad, d)) : indices(itp.coefs, d)
end

function interpolate{TWeights,TC,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(::Type{TWeights}, ::Type{TC}, A, it::IT, gt::GT)
    Apad, Pad = prefilter(TWeights, TC, A, IT, GT)
    BSplineInterpolation(TWeights, Apad, it, gt, Pad)
end
function interpolate{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray, it::IT, gt::GT)
    interpolate(tweight(A), tcoef(A), A, it, gt)
end

# We can't just return a tuple-of-types due to julia #12500
tweight(A::AbstractArray) = Float64
tweight(A::AbstractArray{Float32}) = Float32
tweight(A::AbstractArray{Rational{Int}}) = Rational{Int}
tweight{T<:Integer}(A::AbstractArray{T}) = typeof(float(zero(T)))

tcoef(A::AbstractArray) = eltype(A)
tcoef(A::AbstractArray{Float32}) = Float32
tcoef(A::AbstractArray{Rational{Int}}) = Rational{Int}
tcoef{T<:Integer}(A::AbstractArray{T}) = typeof(float(zero(T)))

interpolate!{TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(::Type{TWeights}, A, it::IT, gt::GT) = BSplineInterpolation(TWeights, prefilter!(TWeights, A, IT, GT), it, gt, Val{0}())
function interpolate!{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray, it::IT, gt::GT)
    interpolate!(tweight(A), A, it, gt)
end

offsetsym(off, d) = off == -1 ? Symbol("ixm_", d) :
                    off ==  0 ? Symbol("ix_", d) :
                    off ==  1 ? Symbol("ixp_", d) :
                    off ==  2 ? Symbol("ixpp_", d) : error("offset $off not recognized")

# Ideally we might want to shift the indices symmetrically, but this
# would introduce an inconsistency, so we just append on the right
@inline indices_removepad(inds::Base.OneTo, pad) = Base.OneTo(length(inds) - 2*pad)
@inline indices_removepad(inds, pad) = oftype(inds, first(inds):last(inds) - 2*pad)
@inline indices_addpad(inds::Base.OneTo, pad) = Base.OneTo(length(inds) + 2*pad)
@inline indices_addpad(inds, pad) = oftype(inds, first(inds):last(inds) + 2*pad)
@inline indices_interior(inds, pad) = first(inds)+pad:last(inds)-pad

"""
    map_repeat(f, a, b)

Equivalent to `(f(a[1], b[1]), f(a[2], b[2]), ...)` if `a` and `b` are
tuples of the same lengths, or `(f(a[1], b), f(a[2], b), ...)` if `b`
is a scalar.
"""
@generated function map_repeat{N}(f, a::NTuple{N,Any}, b::NTuple{N,Any})
    ex = [:(f(a[$i], b[$i])) for i = 1:N]
    quote
        $(Expr(:meta, :inline))
        ($(ex...),)
    end
end
@generated function map_repeat{N}(f, a::NTuple{N,Any}, b)
    ex = [:(f(a[$i], b)) for i = 1:N]
    quote
        $(Expr(:meta, :inline))
        ($(ex...),)
    end
end

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("cubic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")
