export Gridded

struct Gridded{D<:Degree} <: InterpolationType end
Gridded(::D) where {D<:Degree} = Gridded{D}()

griddedtype(::Type{Gridded{D}}) where {D<:Degree} = D

const GridIndex{T} = Union{AbstractVector{T}, Tuple}

# Because Ranges check bounds on getindex, it's actually faster to convert the
# knots to Vectors. It's also good to take a copy, so it doesn't get modified later.
struct GriddedInterpolation{T,N,TCoefs,IT<:DimSpec{Gridded},K<:Tuple{Vararg{Vector}},pad} <: AbstractInterpolation{T,N,IT,OnGrid}
    knots::K
    coefs::Array{TCoefs,N}
end
function GriddedInterpolation(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{TCoefs,N}, ::IT, ::Val{pad}) where {N,TCoefs,TWeights<:Real,IT<:DimSpec{Gridded},pad}
    isleaftype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isleaftype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    knts = mapcollect(knots...)
    for (d,k) in enumerate(knts)
        length(k) == size(A, d) || throw(DimensionMismatch("knot vectors must have the same number of elements as the corresponding dimension of the array"))
        length(k) == 1 && error("dimensions of length 1 not yet supported")  # FIXME
        issorted(k) || error("knot-vectors must be sorted in increasing order")
        iextract(IT, d) != NoInterp || k == collect(1:size(A, d)) || error("knot-vector should be the range 1:$(size(A,d)) for the method Gridded{NoInterp}")
    end
    c = zero(TWeights)
    for _ in 2:N
        c *= c
    end
    T = Base.promote_op(*, typeof(c), TCoefs)

    GriddedInterpolation{T,N,TCoefs,IT,typeof(knts),pad}(knts, A)
end

Base.parent(A::GriddedInterpolation) = A.coefs

# A type-stable version of map(collect, knots)
mapcollect() = ()
@inline mapcollect(k::AbstractVector) = (collect(k),)
@inline mapcollect(k1::AbstractVector, k2::AbstractVector...) = (collect(k1), mapcollect(k2...)...)

# Utilities for working either with scalars or tuples/tuple-types
iextract(::Type{T}, d) where {T<:Gridded} = T
iextract(::Type{T}, d) where {T<:GridType} = T

@generated function size(itp::GriddedInterpolation{T,N,TCoefs,IT,K,pad}, d) where {T,N,TCoefs,IT,K,pad}
    quote
        d <= $N ? size(itp.coefs, d) - 2*padextract($pad, d) : 1
    end
end

function interpolate(::Type{TWeights}, ::Type{TCoefs}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {TWeights,TCoefs,Tel,N,IT<:DimSpec{Gridded}}
    GriddedInterpolation(TWeights, knots, A, it, Val{0}())
end
function interpolate(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {Tel,N,IT<:DimSpec{Gridded}}
    interpolate(tweight(A), tcoef(A), knots, A, it)
end

interpolate!(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {TWeights,Tel,N,IT<:DimSpec{Gridded}} = GriddedInterpolation(TWeights, knots, A, it, Val{0}())
function interpolate!(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {Tel,N,IT<:DimSpec{Gridded}}
    interpolate!(tweight(A), tcoef(A), knots, A, it)
end

lbound(itp::GriddedInterpolation, d) = itp.knots[d][1]
ubound(itp::GriddedInterpolation, d) = itp.knots[d][end]
lbound(itp::GriddedInterpolation, d, inds) = itp.knots[d][1]
ubound(itp::GriddedInterpolation, d, inds) = itp.knots[d][end]

include("constant.jl")
include("linear.jl")
include("indexing.jl")
