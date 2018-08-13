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

struct BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Axs<:Tuple{Vararg{AbstractUnitRange,N}}} <: AbstractInterpolation{T,N,IT,GT}
    coefs::TCoefs
    parentaxes::Axs
    it::IT
    gt::GT
end
function BSplineInterpolation(::Type{TWeights}, A::AbstractArray{Tel,N}, it::IT, gt::GT, axs) where {N,Tel,TWeights<:Real,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    # String interpolation causes allocation, noinline avoids that unless they get called
    @noinline err_concrete(IT) = error("The b-spline type must be a concrete type (was $IT)")
    @noinline warn_concrete(A) = @warn("For performance reasons, consider using an array of a concrete type (typeof(A) == $(typeof(A)))")

    isconcretetype(IT) || err_concrete(IT)
    isconcretetype(typeof(A)) || warn_concrete(A)

    # Compute the output element type when positions have type TWeights
    if isempty(A)
        T = Base.promote_op(*, TWeights, eltype(A))
    else
        T = typeof(zero(TWeights) * first(A))
    end
    BSplineInterpolation{T,N,typeof(A),IT,GT,typeof(axs)}(A, fix_axis.(axs), it, gt)
end

# Utilities for working either with scalars or tuples/tuple-types
iextract(::Type{T}, d) where {T<:BSpline} = T
iextract(t, d) = t.parameters[d]
padding(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}) where {T,N,TCoefs,IT,GT,pad} = pad
padding(itp::AbstractInterpolation) = padding(typeof(itp))
padextract(pad::Integer, d) = pad
padextract(pad::Tuple{Vararg{Integer}}, d) = pad[d]

lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d::Integer) where {T,N,TCoefs,IT} =
    first(axes(itp, d))
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d::Integer) where {T,N,TCoefs,IT} =
    last(axes(itp, d))
lbound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d::Integer) where {T,N,TCoefs,IT} =
    first(axes(itp, d)) - 0.5
ubound(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d::Integer) where {T,N,TCoefs,IT} =
    last(axes(itp, d))+0.5

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

@inline axes(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}) where {T,N,TCoefs,IT,GT,pad} =
    indices_removepad.(axes(itp.coefs), pad)

fix_axis(r::Base.OneTo) = r
fix_axis(r::Base.Slice) = r
fix_axis(r::UnitRange) = Base.Slice(r)
fix_axis(r::AbstractUnitRange) = fix_axis(UnitRange(r))
function axes(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d) where {T,N,TCoefs,IT,GT,pad}
    d <= N ? indices_removepad(axes(itp.coefs, d), padextract(pad, d)) : axes(itp.coefs, d)
end

function interpolate(::Type{TWeights}, ::Type{TC}, A, it::IT, gt::GT) where {TWeights,TC,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    Apad = prefilter(TWeights, TC, A, it, gt)
    BSplineInterpolation(TWeights, Apad, it, gt, axes(A))
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

function interpolate!(::Type{TWeights}, A, it::IT, gt::GT) where {TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    # Set the bounds of the interpolant inward, if necessary
    axsA = axes(A)
    axspad = padded_axes(axsA, it)
    BSplineInterpolation(TWeights, prefilter!(TWeights, A, it, gt), it, gt, fix_axis.(padinset.(axsA, axspad)))
end
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

@static if VERSION < v"0.7.0-DEV.3449"
    lut!(dl, d, du) = lufact!(Tridiagonal(dl, d, du), Val{false})
else
    lut!(dl, d, du) = lu!(Tridiagonal(dl, d, du), Val(false))
end

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("cubic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")

Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT<:Union{BSpline{Linear},BSpline{Constant}}} = A.coefs
Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT} = throw(ArgumentError("The given BSplineInterpolation does not serve as a \"view\" for a parent array. This would only be true for Constant and Linear b-splines."))
