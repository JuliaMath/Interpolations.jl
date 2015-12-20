export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic,
    Cubic

abstract Degree{N} <: Flag

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
iextract{T<:GridType}(::Type{T}, d) = T
iextract(t, d) = t.parameters[d]
padding{T,N,TCoefs,IT,GT,pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}) = pad
padding(itp::AbstractInterpolation) = padding(typeof(itp))
padextract(pad::Integer, d) = pad
padextract(pad::Tuple{Vararg{Integer}}, d) = pad[d]

lbound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d) = 1
ubound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnGrid}, d) = size(itp, d)
lbound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d) = 0.5
ubound{T,N,TCoefs,IT}(itp::BSplineInterpolation{T,N,TCoefs,IT,OnCell}, d) = size(itp, d)+0.5

count_interp_dims{T,N,TCoefs,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType},pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,pad}}, n) = count_interp_dims(IT, n)

function size{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d)
    d <= N ? size(itp.coefs, d) - 2*padextract(pad, d) : 1
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

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("cubic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")
