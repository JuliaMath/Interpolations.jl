module Interpolations

using Base.Cartesian
using WoodburyMatrices

import Base:
    eltype,
    gradient,
    getindex,
    ndims,
    size

export
    Interpolation,
    Interpolation!,
    Constant,
    Linear,
    Quadratic,
    ExtrapError,
    ExtrapNaN,
    ExtrapConstant,
    ExtrapLinear,
    ExtrapReflect,
    ExtrapPeriodic,
    OnCell,
    OnGrid,
    Flat,
    Line,
    Free,
    Periodic,
    Reflect,
    Interior,

    degree,
    boundarycondition,
    gridrepresentation,

    gradient,
    gradient!

abstract Degree{N}

abstract GridRepresentation
immutable OnGrid <: GridRepresentation end
immutable OnCell <: GridRepresentation end

abstract BoundaryCondition
immutable None <: BoundaryCondition end
immutable Flat <: BoundaryCondition end
immutable Line <: BoundaryCondition end
typealias Natural Line
immutable Free <: BoundaryCondition end
immutable Periodic <: BoundaryCondition end
immutable Reflect <: BoundaryCondition end
immutable Interior <: BoundaryCondition end

abstract InterpolationType{D<:Degree,BC<:BoundaryCondition,GR<:GridRepresentation}

include("extrapolation.jl")

abstract AbstractInterpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}
type Interpolation{T,N,TCoefs,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractInterpolation{T,N,IT,EB}
    coefs::Array{TCoefs,N}
end
function Interpolation{N,TCoefs,TWeights<:Real,IT<:InterpolationType,EB<:ExtrapolationBehavior}(::Type{TWeights}, A::AbstractArray{TCoefs,N}, it::IT, ::EB)
    isleaftype(IT) || error("The interpolation type must be a leaf type (was $IT)")

    isleaftype(TCoefs) || warn("For performance reasons, consider using an array of a concrete type T (eltype(A) == $(eltype(A)))")

    c = one(TWeights)
    for _ in 2:N
        c *= c
    end
    T = typeof(c * one(TCoefs))

    Interpolation{T,N,TCoefs,IT,EB}(prefilter(TWeights,A,it))
end
Interpolation(A::AbstractArray, it::InterpolationType, eb::ExtrapolationBehavior) = Interpolation(Float64, A, it, eb)
Interpolation(A::AbstractArray{Float32}, it::InterpolationType, eb::ExtrapolationBehavior) = Interpolation(Float32, A, it, eb)
Interpolation(A::AbstractArray{Rational{Int}}, it::InterpolationType, eb::ExtrapolationBehavior) = Interpolation(Rational{Int}, A, it, eb)

# In-place (destructive) interpolation
function Interpolation!{T<:FloatingPoint,N,D<:Degree,G<:GridRepresentation,EB<:ExtrapolationBehavior}(A::AbstractArray{T,N}, it::InterpolationType{D,Interior,G}, ::EB)
    IT = typeof(it)
    isleaftype(typeof(IT)) || error("The interpolation type must be a leaf type (was $IT)")
    isleaftype(T) || error("Must be an array of a concrete type T (eltype(A) == $(eltype(A)))")
    Interpolation{T,N,T,IT,EB}(prefilter!(T,A,0,it))
end

# Unless otherwise specified, use coefficients as they are, i.e. without prefiltering
# However, all prefilters copy the array, so do that here as well
prefilter{TWeights,T,N,IT<:InterpolationType}(::Type{TWeights}, A::AbstractArray{T,N}, ::IT) = copy(A)

size{T,N}(itp::Interpolation{T,N}, d::Integer) = d > N ? 1 : size(itp.coefs, d) - 2*padding(interpolationtype(itp))
size(itp::AbstractInterpolation) = tuple([size(itp,i) for i in 1:ndims(itp)]...)
ndims(itp::Interpolation) = ndims(itp.coefs)

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")

interpolationtype{T,N,TCoefs,IT<:InterpolationType}(itp::Interpolation{T,N,TCoefs,IT}) = IT()

boundarycondition{D,BC<:BoundaryCondition}(::InterpolationType{D,BC}) = BC()
gridrepresentation{D,BC,GR<:GridRepresentation}(::InterpolationType{D,BC,GR}) = GR()
degree{D<:Degree,BC,GR}(::InterpolationType{D,BC,GR}) = D()

# If nothing else is specified, don't pad at all
padding(::InterpolationType) = 0

for op in (:boundarycondition, :gridrepresentation, :degree, :padding)
    eval(:($(op)(itp::Interpolation) = $(op)(interpolationtype(itp))))
end

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("filter1d.jl")

function padded_index(sz::Tuple, pad)
    sz = Int[s+2pad for s in sz]
    N = length(sz)
    ind = cell(N)
    for i in 1:N
        ind[i] = 1+pad:sz[i]-pad
    end
    ind,sz
end
function copy_with_padding(A, it::InterpolationType)
    pad = padding(it)
    ind,sz = padded_index(size(A), pad)
    coefs = zeros(eltype(A), sz...)
    coefs[ind...] = A
    coefs, pad
end

supported(::Type{Constant{OnGrid}}) = true
supported(::Type{Constant{OnCell}}) = true
supported(::Type{Linear{OnGrid}}) = true
supported(::Type{Linear{OnCell}}) = true
supported(::Type{Quadratic{Flat,OnCell}}) = true
supported(::Type{Quadratic{Flat,OnGrid}}) = true
supported(::Type{Quadratic{Reflect,OnCell}}) = true
supported(::Type{Quadratic{Reflect,OnGrid}}) = true
supported(::Type{Quadratic{Line,OnGrid}}) = true
supported(::Type{Quadratic{Line,OnCell}}) = true
supported(::Type{Quadratic{Free,OnGrid}}) = true
supported(::Type{Quadratic{Free,OnCell}}) = true
supported(::Type{Quadratic{Periodic,OnGrid}}) = true
supported(::Type{Quadratic{Periodic,OnCell}}) = true
supported(::Any) = false

function getindex_impl(N, IT::Type, EB::Type)
    it = IT()
    eb = EB()
    gr = gridrepresentation(it)
    quote
        @nexprs $N d->(x_d = x[d])

        # Handle extrapolation, by either throwing an error,
        # returning a value, or shifting x to somewhere inside
        # the domain
        $(extrap_transform_x(gr,eb,N,it))

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(it, N))

        # Calculate coefficient weights based on fx
        $(coefficients(it, N))

        # Generate the indexing expression
        @inbounds ret = $(index_gen(degree(it), N))
        ret
    end
end

# Resolve ambiguity with linear indexing
getindex{T,N,TCoefs,IT,EB}(itp::Interpolation{T,N,TCoefs,IT,EB}, x::Real) =
    error("Linear indexing is not supported for interpolation objects")
@generated function getindex{T,TCoefs,IT,EB}(itp::Interpolation{T,1,TCoefs,IT,EB}, x::Real)
    getindex_impl(1, IT, EB)
end

# Resolve ambiguity with real multilinear indexing
@generated function getindex{T,N,TCoefs,IT,EB}(itp::Interpolation{T,N,TCoefs,IT,EB}, x::Real...)
    n = length(x)
    n == N || return :(error("Must index ", $N, "-dimensional interpolation objects with ", $(nindexes(N))))
    getindex_impl(N, IT, EB)
end

# Resolve ambiguity with general array indexing,
#   getindex{T,N}(A::AbstractArray{T,N}, I::AbstractArray{T,N})
function getindex{T,N,TCoefs,IT,EB}(itp::Interpolation{T,N,TCoefs,IT,EB}, I::AbstractArray{T})
    error("Array indexing is not defined for interpolation objects.")
end

# Resolve ambiguity with colon indexing for 1D interpolations
#   getindex{T}(A::AbstractArray{T,1}, C::Colon)
function getindex{T,TCoefs,IT,EB}(itp::Interpolation{T,1,TCoefs,IT,EB}, c::Colon)
    error("Colon indexing is not supported for interpolation objects")
end

# Allow multilinear indexing with any types
@generated function getindex{T,N,TCoefs,TIndex,IT,EB}(itp::Interpolation{T,N,TCoefs,IT,EB}, x::TIndex...)
    n = length(x)
    n == N || return :(error("Must index ", $N, "-dimensional interpolation objects with ", $(nindexes(N))))
    getindex_impl(N, IT, EB)
end

# Define in-place gradient calculation
#   gradient!(g::Vector, itp::Interpolation, NTuple{N}...)
@generated function gradient!{T,N,TCoefs,TOut,IT,EB}(g::Vector{TOut}, itp::Interpolation{T,N,TCoefs,IT,EB}, x...)
    n = length(x)
    n == N || return :(error("Must index ", $N, "-dimensional interpolation objects with ", $(nindexes(N))))
    it = IT()
    eb = EB()
    gr = gridrepresentation(it)
    quote
        length(g) == $N || error("g must be an array with exactly N elements (length(g) == "*string(length(g))*", N == "*string($N)*")")
        @nexprs $N d->(x_d = x[d])
        $(extrap_transform_x(gr,eb,N,it))
        $(define_indices(it,N))
        @nexprs $N dim->begin
            @nexprs $N d->begin
                (d==dim ? $(gradient_coefficients(it,N,:d))
                        : $(coefficients(it,N,:d)))
            end

            @inbounds g[dim] = $(index_gen(degree(it),N))
        end
        g
    end
end

gradient{T,N}(itp::Interpolation{T,N}, x...) = gradient!(Array(T,N), itp, x...)

function prefilter{TWeights,TCoefs,N,IT<:Quadratic}(::Type{TWeights}, A::Array{TCoefs,N}, it::IT)
    ret, pad = copy_with_padding(A, it)
    prefilter!(TWeights, ret, pad, it)
end

function prefilter!{TWeights,TCoefs,N,IT<:Quadratic}(::Type{TWeights}, ret::Array{TCoefs,N}, pad::Integer, it::IT)
    sz = size(ret)
    first = true
    for dim in 1:N
        M, b = prefiltering_system(TWeights, TCoefs, sz[dim], it)
        if !isa(M, Woodbury)
            A_ldiv_B_md!(ret, M, ret, dim, b)
        else
            if first
                buf = Array(eltype(ret), length(ret,))
                shape = sz
                retrs = reshape(ret, shape)  # for type-stability against a future julia #10507
                first = false
            end
            bufrs = reshape(buf, shape)
            filter_dim1!(bufrs, M, retrs, b)
            shape = (sz[dim+1:end]..., sz[1:dim]...)
            retrs = reshape(ret, shape)
            permutedims!(retrs, bufrs, ((2:N)..., 1))
        end
    end
    ret
end

nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"

end # module
