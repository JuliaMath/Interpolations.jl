module Interpolations

using Base.Cartesian
using Compat

import Base: size, eltype, getindex, ndims

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

    degree,
    boundarycondition,
    gridrepresentation

abstract Degree{N}

abstract GridRepresentation
type OnGrid <: GridRepresentation end
type OnCell <: GridRepresentation end

abstract BoundaryCondition
type None <: BoundaryCondition end
type Flat <: BoundaryCondition end
type Line <: BoundaryCondition end
typealias Natural Line
type Free <: BoundaryCondition end
type Periodic <: BoundaryCondition end

abstract InterpolationType{D<:Degree,BC<:BoundaryCondition,GR<:GridRepresentation}

include("extrapolation.jl")

abstract AbstractInterpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}
type Interpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractInterpolation{T,N,IT,EB}
    coefs::Array{T,N}
end
function Interpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior}(A::Array{T,N}, it::IT, ::EB)
    isleaftype(IT) || error("The interpolation type must be a leaf type (was $IT)")
    
    isleaftype(T) || warn("For performance reasons, consider using an array of a concrete type T (eltype(A) == $(eltype(A)))")

    Interpolation{T,N,IT,EB}(prefilter(A,it))
end

# Unless otherwise specified, use coefficients as they are, i.e. without prefiltering
# However, all prefilters copy the array, so do that here as well
prefilter{T,N,IT<:InterpolationType}(A::AbstractArray{T,N}, ::IT) = copy(A)

size{T,N,IT<:InterpolationType}(itp::Interpolation{T,N,IT}, d::Integer) =
    size(itp.coefs, d) - 2*padding(IT())
size(itp::Interpolation) = tuple([size(itp,i) for i in 1:ndims(itp)]...)
ndims(itp::Interpolation) = ndims(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

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
padding(::InterpolationType) = 0
padding{T,N,IT<:InterpolationType}(::Interpolation{T,N,IT}) = padding(IT())

function pad_size_and_index(sz::Tuple, pad)
    sz = Int[s+2pad for s in sz]
    N = length(sz)
    ind = cell(N)
    for i in 1:N
        ind[i] = 1+pad:sz[i]-pad
    end
    sz, ind
end
function copy_with_padding(A, it::InterpolationType)
    pad = padding(it)
    sz,ind = pad_size_and_index(size(A), pad)
    coefs = fill(convert(eltype(A), 0), sz...)
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
            :(promote_type(T, x...)),
            :(getindex{T,N}(itp::Interpolation{T,N,$IT,$EB}, x::NTuple{N,Real}...)), 
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
    end
end

# This creates prefilter specializations for all interpolation types that need them
for IT in (
        Quadratic{Flat,OnCell},
        Quadratic{Flat,OnGrid},
        Quadratic{Line,OnGrid},
        Quadratic{Line,OnCell},
        Quadratic{Free,OnGrid},
        Quadratic{Free,OnCell},
        Quadratic{Periodic,OnGrid},
        Quadratic{Periodic,OnCell},
    )
    @ngenerate N promote_type_grid(T, x...) function prefilter{T,N}(A::Array{T,N},it::IT)
        ret, pad = copy_with_padding(A,it)

        szs = collect(size(ret))
        strds = collect(strides(ret))

        for dim in 1:N
            n = szs[dim]
            szs[dim] = 1

            M, b = prefiltering_system(eltype(A), n, it)

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
