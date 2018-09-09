export Gridded

struct Gridded{D<:Degree} <: InterpolationType
    degree::D
end

function Base.show(io::IO, g::Gridded)
    print(io, "Gridded(")
    show(io, degree(g))
    print(io, ')')
end

const GridIndex{T} = Union{AbstractVector{T}, Tuple}

struct GriddedInterpolation{T,N,TCoefs,IT<:DimSpec{Gridded},K<:Tuple{Vararg{AbstractVector}}} <: AbstractInterpolation{T,N,IT}
    knots::K
    coefs::Array{TCoefs,N}
    it::IT
end
function GriddedInterpolation(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{TCoefs,N}, it::IT) where {N,TCoefs,TWeights<:Real,IT<:DimSpec{Gridded},pad}
    isconcretetype(IT) || error("The b-spline type must be a leaf type (was $IT)")
    isconcretetype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    check_gridded(it, knots, axes(A))
    c = zero(TWeights)
    if isempty(A)
        T = Base.promote_op(*, typeof(c), eltype(A))
    else
        T = typeof(c * first(A))
    end
    GriddedInterpolation{T,N,TCoefs,IT,typeof(knots)}(knots, A, it)
end

@inline function check_gridded(itpflag, knots, axs)
    flag, ax1, k1 = getfirst(itpflag), axs[1], knots[1]
    if flag isa NoInterp
        k1 == ax1 || error("for NoInterp knot vector should be $ax1, got $k1")
    else
        axes(k1, 1) == ax1 || throw(DimensionMismatch("knot vectors must have the same axes as the corresponding dimension of the array"))
    end
    degree(flag) isa Union{NoInterp,Constant,Linear} || error("only Linear, Constant, and NoInterp supported, got $flag")
    length(k1) == 1 && error("dimensions of length 1 not yet supported")  # FIXME
    issorted(k1) || error("knot-vectors must be sorted in increasing order")
    check_gridded(getrest(itpflag), Base.tail(knots), Base.tail(axs))
end
check_gridded(::Any, ::Tuple{}, ::Tuple{}) = nothing
degree(flag::Gridded) = flag.degree

Base.parent(A::GriddedInterpolation) = A.coefs
coefficients(A::GriddedInterpolation) = A.coefs

size(A::GriddedInterpolation) = size(A.coefs)
axes(A::GriddedInterpolation) = axes(A.coefs)

itpflag(A::GriddedInterpolation) = A.it

function interpolate(::Type{TWeights}, ::Type{TCoefs}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {TWeights,TCoefs,Tel,N,IT<:DimSpec{Gridded}}
    GriddedInterpolation(TWeights, knots, A, it)
end
function interpolate(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {Tel,N,IT<:DimSpec{Gridded}}
    interpolate(tweight(A), tcoef(A), knots, A, it)
end

interpolate!(::Type{TWeights}, knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {TWeights,Tel,N,IT<:DimSpec{Gridded}} =
    GriddedInterpolation(TWeights, knots, A, it)
function interpolate!(knots::NTuple{N,GridIndex}, A::AbstractArray{Tel,N}, it::IT) where {Tel,N,IT<:DimSpec{Gridded}}
    interpolate!(tweight(A), tcoef(A), knots, A, it)
end

lbounds(itp::GriddedInterpolation) = first.(itp.knots)
ubounds(itp::GriddedInterpolation) = last.(itp.knots)

include("constant.jl")
include("linear.jl")
include("indexing.jl")
