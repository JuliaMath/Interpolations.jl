module Interpolations

using Base.Cartesian
using Compat

import Base:
    eltype,
    gradient,
    getindex,
    ndims,
    size

export 
    Interpolation,
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

abstract InterpolationType{D<:Degree,BC<:BoundaryCondition,GR<:GridRepresentation}

include("extrapolation.jl")

abstract AbstractInterpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}
type Interpolation{TEl,N,TCoefs<:Real,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractInterpolation{TEl,N,IT,EB}
    coefs::Array{TCoefs,N}
end
function Interpolation{TIn,N,TCoefs,IT<:InterpolationType,EB<:ExtrapolationBehavior}(::Type{TCoefs}, A::AbstractArray{TIn,N}, it::IT, ::EB)
    isleaftype(IT) || error("The interpolation type must be a leaf type (was $IT)")
    
    isleaftype(TIn) || warn("For performance reasons, consider using an array of a concrete type T (eltype(A) == $(eltype(A)))")

    c = one(TCoefs)
    for _ in 2:N
        c *= c
    end
    TEl = typeof(c)

    Interpolation{TEl,N,TCoefs,IT,EB}(prefilter(TCoefs,A,it))
end
Interpolation{TIn,N,IT<:InterpolationType,EB<:ExtrapolationBehavior}(A::AbstractArray{TIn,N}, it::IT, eb::EB) = Interpolation(TIn, A, it, eb)

# Unless otherwise specified, use coefficients as they are, i.e. without prefiltering
# However, all prefilters copy the array, so do that here as well
# We also ensure that the coefficient array is of the correct type
prefilter{T,N,TCoefs,IT<:InterpolationType}(::Type{TCoefs}, A::AbstractArray{T,N}, ::IT) = copy!(Array(TCoefs,size(A)...), A)

size{T,N,TCoefs,IT<:InterpolationType}(itp::Interpolation{T,N,TCoefs,IT}, d::Integer) = size(itp.coefs, d) - 2*padding(TCoefs,IT())
size(itp::AbstractInterpolation) = tuple([size(itp,i) for i in 1:ndims(itp)]...)
ndims(itp::Interpolation) = ndims(itp.coefs)
eltype{T}(itp::Interpolation{T}) = T
coeftype(itp::Interpolation) = eltype(itp.coefs)

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")

boundarycondition{D,BC<:BoundaryCondition}(::InterpolationType{D,BC}) = BC()
boundarycondition{T,N,IT}(::Interpolation{T,N,IT}) = boundarycondition(IT())
gridrepresentation{D,BC,GR<:GridRepresentation}(::InterpolationType{D,BC,GR}) = GR()
gridrepresentation{T,N,IT}(::Interpolation{T,N,IT}) = gridrepresentation(IT())
degree{D<:Degree,BC,GR}(::InterpolationType{D,BC,GR}) = D()
degree{T,N,IT}(::Interpolation{T,N,IT}) = degree(IT())

include("constant.jl")
include("linear.jl")
include("quadratic.jl")

# If nothing else is specified, don't pad at all
padding{T}(::Type{T}, ::InterpolationType) = zero(T)
padding{T,N,IT<:InterpolationType}(::Interpolation{T,N,IT}) = padding(T,IT())

function pad_size_and_index(sz::Tuple, pad)
    sz = Int[s+2pad for s in sz]
    N = length(sz)
    ind = cell(N)
    for i in 1:N
        ind[i] = 1+pad:sz[i]-pad
    end
    sz, ind
end
function copy_with_padding(TCoefs, A, it::InterpolationType)
    pad = padding(it)
    sz,ind = pad_size_and_index(size(A), pad)
    coefs = fill!(Array(TCoefs, sz...), 0)
    coefs[ind...] = A
    coefs, pad
end

# This creates getindex methods for all supported combinations
for IT in (
        Constant{OnGrid},
        Constant{OnCell},
        Linear{OnGrid},
        Linear{OnCell},
        Quadratic{Flat,OnCell},
        Quadratic{Flat,OnGrid},
        Quadratic{Reflect,OnCell},
        Quadratic{Reflect,OnGrid},
        Quadratic{Line,OnGrid},
        Quadratic{Line,OnCell},
        Quadratic{Free,OnGrid},
        Quadratic{Free,OnCell},
        Quadratic{Periodic,OnGrid},
        Quadratic{Periodic,OnCell},
    )
    for EB in (
            ExtrapError,
            ExtrapNaN,
            ExtrapConstant,
            ExtrapLinear,
            ExtrapReflect,
            ExtrapPeriodic,
        )
        it = IT()
        eb = EB()
        gr = gridrepresentation(it)
        eval(ngenerate(
            :N,
            :T,
            :(getindex{T,N,TCoefs<:Real}(itp::Interpolation{T,N,TCoefs,$IT,$EB}, x::NTuple{N,TCoefs}...)), 
            N->quote
                # Handle extrapolation, by either throwing an error,
                # returning a value, or shifting x to somewhere inside
                # the domain
                $(extrap_transform_x(gr,eb,N))

                # Calculate the indices of all coefficients that will be used
                # and define fx = x - xi in each dimension
                $(define_indices(it, N))

                # Calculate coefficient weights based on fx
                $(coefficients(it, N))

                # Generate the indexing expression
                @inbounds ret = $(index_gen(degree(it), N))
                ret
            end
        ))
        eval(ngenerate(:N, :T, :(getindex{T,N,TCoefs<:Real}(itp::Interpolation{T,N,TCoefs,$IT,$EB}, xs::Real...)),
            N->:(getindex(itp, [convert(TCoefs, x) for x in xs]...))
        ))

        eval(ngenerate(
            :N,
            :T,
            :(gradient!{T,N,TCoefs<:Real}(g::Array{T,1}, itp::Interpolation{T,N,TCoefs,$IT,$EB}, x::NTuple{N,TCoefs}...)),
            N->quote
                $(extrap_transform_x(gr,eb,N))
                $(define_indices(it,N))
                @nexprs $N dim->begin
                    @nexprs $N d->begin
                        (d==dim
                            ? $(gradient_coefficients(it,N,:d))
                            : $(coefficients(it,N,:d)))
                    end

                    @inbounds g[dim] = $(index_gen(degree(it),N))
                end
                g
            end
        ))
        eval(ngenerate(:N, :T, :(gradient!{T,N,TCoefs<:Real}(g::Array{T,1}, itp::Interpolation{T,N,TCoefs,$IT,$EB}, xs::Real...)),
            N->:(gradient!(g, itp, [convert(TCoefs, x) for x in xs]...))
        ))
    end
end

gradient{T,N}(itp::Interpolation{T,N}, x...) = gradient!(Array(T,N),itp,x...)
gradient1{T}(itp::Interpolation{T,1}, x) = gradient!(Array(T,1),itp,x)[1]

# This creates prefilter specializations for all interpolation types that need them
for IT in (
        Quadratic{Flat,OnCell},
        Quadratic{Flat,OnGrid},
        Quadratic{Reflect,OnCell},
        Quadratic{Reflect,OnGrid},
        Quadratic{Line,OnGrid},
        Quadratic{Line,OnCell},
        Quadratic{Free,OnGrid},
        Quadratic{Free,OnCell},
        Quadratic{Periodic,OnGrid},
        Quadratic{Periodic,OnCell},
    )
    @ngenerate N Array{TCoefs,N} function prefilter{T,N,TCoefs}(::Type{TCoefs}, A::Array{T,N},it::IT)
        ret, pad = copy_with_padding(TCoefs, A,it)

        szs = collect(size(ret))
        strds = collect(strides(ret))

        for dim in 1:N
            n = szs[dim]
            szs[dim] = 1

            M, b = prefiltering_system(TCoefs, n, it)

            @nloops N i d->1:szs[d] begin
                cc = @ntuple N i

                strt_diff = sum([(cc[i]-1)*strds[i] for i in 1:length(cc)])
                strt = 1 + strt_diff
                rng = range(strt, strds[dim], n)

                bdiff = ret[rng]
                b += bdiff
                ret[rng] = M \ b
                b -= bdiff
            end
            szs[dim] = n
        end
        ret
    end
end

end # module
