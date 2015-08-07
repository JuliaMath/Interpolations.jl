export Gridded

immutable Gridded{D<:Degree} <: InterpolationType end
Gridded{D<:Degree}(::Type{D}) = Gridded{D}

griddedtype{D<:Degree}(::Type{Gridded{D}}) = D

typealias GridIndex{T} Union(AbstractVector{T}, Tuple)

immutable GriddedInterpolation{T,N,TCoefs,IT<:DimSpec{Gridded},K<:Tuple{Vararg{GridIndex}},pad} <: AbstractInterpolation{T,N,IT,OnGrid}
    knots::K
    coefs::Array{TCoefs,N}
end
function GriddedInterpolation{N,TCoefs,TWeights<:Real,IT<:DimSpec{Gridded},pad}(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{TCoefs,N}, ::Type{IT}, ::Val{pad})
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    for (d,k) in enumerate(knots)
        length(k) == size(A, d) || throw(DimensionMismatch("knot vectors must have the same number of elements as the corresponding dimension of the array"))
        length(k) == 1 && error("dimensions of length 1 not yet supported")  # FIXME
        issorted(k) || error("knot-vectors must be sorted in increasing order")
    end
    c = one(TWeights)
    for _ in 2:N
        c *= c
    end
    T = typeof(c * one(TCoefs))

    GriddedInterpolation{T,N,TCoefs,IT,typeof(knots),pad}(knots, A)
end

# Utilities for working either with scalars or tuples/tuple-types
iextract{T<:Gridded}(::Type{T}, d) = T
iextract{T<:GridType}(::Type{T}, d) = T

@generated function size{T,N,TCoefs,IT,K,pad}(itp::GriddedInterpolation{T,N,TCoefs,IT,K,pad}, d)
    quote
        d <= $N ? size(itp.coefs, d) - 2*padextract($pad, d) : 1
    end
end

function interpolate{TWeights,TCoefs,Tel,N,IT<:DimSpec{Gridded}}(::Type{TWeights}, ::Type{TCoefs}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, ::Type{IT})
    GriddedInterpolation(TWeights, knots, A, IT, Val{0}())
end
function interpolate{Tel,N,IT<:DimSpec{Gridded}}(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, ::Type{IT})
    interpolate(tweight(A), tcoef(A), knots, A, IT)
end

interpolate!{TWeights,Tel,N,IT<:DimSpec{Gridded}}(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, ::Type{IT}) = GriddedInterpolation(TWeights, knots, A, IT, Val{0}())
function interpolate!{Tel,N,IT<:DimSpec{Gridded}}(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, ::Type{IT})
    interpolate!(tweight(A), tcoef(A), knots, A, IT)
end

include("nointerp.jl")
include("constant.jl")
include("linear.jl")
include("indexing.jl")
