export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic

abstract Degree{N}

immutable BSpline{D<:Degree} <: InterpolationType end
BSpline{D<:Degree}(::Type{D}) = BSpline{D}

bsplinetype{D<:Degree}(::Type{BSpline{D}}) = D

immutable BSplineInterpolation{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad} <: AbstractInterpolation{T,N,IT,GT}
    coefs::Array{TCoefs,N}
end
function BSplineInterpolation{N,TCoefs,TWeights<:Real,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},pad}(::Type{TWeights}, A::AbstractArray{TCoefs,N}, ::Type{IT}, ::Type{GT}, ::Val{pad})
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    c = one(TWeights)
    for _ in 2:N
        c *= c
    end
    T = typeof(c * one(TCoefs))

    BSplineInterpolation{T,N,TCoefs,IT,GT,pad}(A)
end

# Utilities for working either with scalars or tuples/tuple-types
iextract{T<:BSpline}(::Type{T}, d) = T
iextract{T<:GridType}(::Type{T}, d) = T
iextract(t, d) = t.parameters[d]
padextract(pad::Integer, d) = pad
padextract(pad::Tuple{Vararg{Integer}}, d) = pad[d]

@generated function size{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, d)
    quote
        d <= $N ? size(itp.coefs, d) - 2*padextract($pad, d) : 1
    end
end

function interpolate{TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT})
    Apad, Pad = prefilter(TWeights, A, IT, GT)
    BSplineInterpolation(TWeights, Apad, IT, GT, Pad)
end
interpolate{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray, ::Type{IT}, ::Type{GT}) = interpolate(Float64, A, IT, GT)
interpolate{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray{Float32}, ::Type{IT}, ::Type{GT}) = interpolate(Float32, A, IT, GT)
interpolate{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray{Rational{Int}}, ::Type{IT}, ::Type{GT}) = interpolate(Rational{Int}, A, IT, GT)

interpolate!{TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = BSplineInterpolation(TWeights, prefilter!(TWeights, A, IT, GT), IT, GT, Val{0}())
interpolate!{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray, ::Type{IT}, ::Type{GT}) = interpolate!(Float64, A, IT, GT)
interpolate!{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray{Float32}, ::Type{IT}, ::Type{GT}) = interpolate!(Float32, A, IT, GT)
interpolate!{IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}(A::AbstractArray{Rational{Int}}, ::Type{IT}, ::Type{GT}) = interpolate!(Rational{Int}, A, IT, GT)

define_indices{IT}(::Type{IT}, N, pad) = Expr(:block, Expr[define_indices_d(iextract(IT, d), d, padextract(pad, d)) for d = 1:N]...)

coefficients{IT}(::Type{IT}, N) = Expr(:block, Expr[coefficients(iextract(IT, d), N, d) for d = 1:N]...)

function gradient_coefficients{IT<:DimSpec{BSpline}}(::Type{IT}, N, dim)
    exs = Expr[d==dim ? gradient_coefficients(iextract(IT, dim), d) :
                        coefficients(iextract(IT, d), N, d) for d = 1:N]
    Expr(:block, exs...)
end

index_gen{IT}(::Type{IT}, N::Integer, offsets...) = index_gen(iextract(IT, min(length(offsets)+1, N)), IT, N, offsets...)

@generated function gradient{T,N,TCoefs,IT,GT,pad}(itp::BSplineInterpolation{T,N,TCoefs,IT,GT,pad}, xs...)
    n = count_interp_dims(IT, N)
    :(gradient!(Array(T,$n), itp, xs...))
end

include("nointerp.jl")
include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")
