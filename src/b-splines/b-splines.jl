export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic,
    Cubic

abstract type Degree{N} <: Flag end

struct BSpline{D<:Degree} <: InterpolationType end
BSpline(::D) where {D<:Degree} = BSpline{D}()

bsplinetype(::Type{BSpline{D}}) where {D<:Degree} = D

struct BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad} <: AbstractInterpolation{T,N,IT,GT}
    coefs::TCoefs
end
function BSplineInterpolation(::Type{TWeights}, A::AbstractArray{Tel,N}, ::IT, ::GT, ::Val{pad}) where {N,Tel,TWeights<:Real,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad}
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(typeof(A)) || warn("For performance reasons, consider using an array of a concrete type (typeof(A) == $(typeof(A)))")

    c = zero(TWeights)
    for _ in 2:N
        c *= c
    end
    if isempty(A)
        T = Base.promote_op(*, typeof(c), eltype(A))
    else
        T = typeof(c * first(A))
    end
    BSplineInterpolation{T,N,typeof(A),IT,GT,pad}(A)
end

# Utilities for working either with scalars or tuples/tuple-types
iextract(::Type{T}, d) where {T<:BSpline} = T
iextract(t, d) = t.parameters[d]
padding(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}) where {T,N,TCoefs,IT,GT,pad} = pad
padding(itp::AbstractInterpolation) = padding(typeof(itp))
padextract(pad::Integer, d) = pad
padextract(pad::Tuple{Vararg{Integer}}, d) = pad[d]

lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d::Integer) where {T,N,TCoefs,IT} =
    first(indices(itp, d))
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d::Integer) where {T,N,TCoefs,IT} =
    last(indices(itp, d))
lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d::Integer) where {T,N,TCoefs,IT} =
    first(indices(itp, d)) - 0.5
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d::Integer) where {T,N,TCoefs,IT} =
    last(indices(itp, d))+0.5

lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d, inds) where {T,N,TCoefs,IT} =
    first(inds)
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d, inds) where {T,N,TCoefs,IT} =
    last(inds)
lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d, inds) where {T,N,TCoefs,IT} =
    first(inds) - 0.5
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d, inds) where {T,N,TCoefs,IT} =
    last(inds)+0.5

count_interp_dims(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}, n) where {T,N,TCoefs,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType},pad} = count_interp_dims(IT, n)

function size(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d) where {T,N,TCoefs,IT,GT,pad}
    d <= N ? size(itp.coefs, d) - 2*padextract(pad, d) : 1
end

@inline indices(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}) where {T,N,TCoefs,IT,GT,pad} =
    indices_removepad.(indices(itp.coefs), pad)

function indices(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d) where {T,N,TCoefs,IT,GT,pad}
    d <= N ? indices_removepad(indices(itp.coefs, d), padextract(pad, d)) : indices(itp.coefs, d)
end

function interpolate(::Type{TWeights}, ::Type{TC}, A, it::IT, gt::GT) where {TWeights,TC,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    Apad, Pad = prefilter(TWeights, TC, A, IT, GT)
    BSplineInterpolation(TWeights, Apad, it, gt, Pad)
end
function interpolate(A::AbstractArray, it::IT, gt::GT) where {IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    interpolate(tweight(A), tcoef(A), A, it, gt)
end

# We can't just return a tuple-of-types due to julia #12500
tweight(A::AbstractArray) = Float64
tweight(A::AbstractArray{Float32}) = Float32
tweight(A::AbstractArray{Rational{Int}}) = Rational{Int}
tweight(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

tcoef(A::AbstractArray) = eltype(A)
tcoef(A::AbstractArray{Float32}) = Float32
tcoef(A::AbstractArray{Rational{Int}}) = Rational{Int}
tcoef(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

interpolate!(::Type{TWeights}, A, it::IT, gt::GT) where {TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}} = BSplineInterpolation(TWeights, prefilter!(TWeights, A, IT, GT), it, gt, Val{0}())
function interpolate!(A::AbstractArray, it::IT, gt::GT) where {IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
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

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("cubic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")

Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT<:Union{BSpline{Linear},BSpline{Constant}}} = A.coefs
Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT} = throw(ArgumentError("The given BSplineInterpolation does not serve as a \"view\" for a parent array. This would only be true for Constant and Linear b-splines."))
