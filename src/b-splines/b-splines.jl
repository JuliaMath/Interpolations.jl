export
    BSpline,

    Constant,
    Linear

abstract Degree{N}

immutable BSpline{D<:Degree} <: InterpolationType end
BSpline{D<:Degree}(::Type{D}) = BSpline{D}

type BSplineInterpolation{T,N,TCoefs,IT<:BSpline,GT<:GridType} <: AbstractInterpolation{T,N,IT,GT}
    coefs::Array{TCoefs,N}
end
function BSplineInterpolation{N,TCoefs,TWeights<:Real,IT<:BSpline,GT<:GridType}(::Type{TWeights}, A::AbstractArray{TCoefs,N}, ::Type{IT}, ::Type{GT})
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    c = one(TWeights)
    for _ in 2:N
        c *= c
    end
    T = typeof(c * one(TCoefs))

    BSplineInterpolation{T,N,TCoefs,IT,GT}(A)
end

size(itp::BSplineInterpolation, d) = size(itp.coefs, d)

interpolate{TWeights,IT<:BSpline,GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = BSplineInterpolation(TWeights, prefilter(TWeights, A, IT, GT), IT, GT)
interpolate{IT<:BSpline,GT<:GridType}(A::AbstractArray, ::Type{IT}, ::Type{GT}) = interpolate(Float64, A, IT, GT)
interpolate{IT<:BSpline,GT<:GridType}(A::AbstractArray{Float32}, ::Type{IT}, ::Type{GT}) = interpolate(Float32, A, IT, GT)
interpolate{IT<:BSpline,GT<:GridType}(A::AbstractArray{Rational{Int}}, ::Type{IT}, ::Type{GT}) = interpolate(Rational{Int}, A, IT, GT)

prefilter{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A

include("constant.jl")
include("linear.jl")
# include("quadratic.jl")
include("indexing.jl")
