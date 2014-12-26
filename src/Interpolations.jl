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
    LinearBC,
    Free,
    Periodic

abstract Degree{N}

abstract GridRepresentation
type OnGrid <: GridRepresentation end
type OnCell <: GridRepresentation end

abstract BoundaryCondition
type None <: BoundaryCondition end
type Flat <: BoundaryCondition end
type LinearBC <: BoundaryCondition end
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
size(itp::Interpolation) = tuple(collect([size(itp,i) for i in 1:ndims(itp)]))
ndims(itp::Interpolation) = ndims(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")

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
function similar_with_padding(A, it::InterpolationType)
    pad = padding(it)
    coefs = Array(eltype(A), [size(A,i)+2pad for i in 1:ndims(A)]...)
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
        Quadratic{LinearBC,OnGrid},
        Quadratic{LinearBC,OnCell},
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
        eval(:(function getindex{T}(itp::Interpolation{T,1,$IT,$EB}, x::Real, d)
            d == 1 || throw(BoundsError())
            itp[x]
        end))

        it = IT()
        eb = EB()
        gr = gridrepresentation(it)
        eval(ngenerate(
            :N,
            :(promote_type(T, x...)),
            :(getindex{T,N}(itp::Interpolation{T,N,$IT,$EB}, x::NTuple{N,Real}...)), 
            N->quote
                $(extrap_gen(gr,eb,N))

                # If the boundary condition mandates separate treatment, this is done
                # by bc_gen.
                # Given an interpolation object itp, with N dimensions, and a coordinate
                # x_d, it should define ix_d such that all the coefficients required by
                # the interpolation will be inbounds.
                $(bc_gen(it, N))

                # indices calculates the indices required for this interpolation degree,
                # based on ix_d defined by bc_gen(), as well as the distance fx_d from 
                # the cell index ix_d to the interpolation coordinate x_d
                $(indices(it, N))

                # These coefficients determine the interpolation basis expansion
                $(coefficients(it, N))

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
        Quadratic{LinearBC,OnGrid},
        Quadratic{LinearBC,OnCell},
        Quadratic{Free,OnGrid},
        Quadratic{Free,OnCell},
        Quadratic{Periodic,OnGrid},
        Quadratic{Periodic,OnCell},
    )
    @ngenerate N promote_type_grid(T, x...) function prefilter{T,N}(A::Array{T,N},it::IT)
        ret, pad = similar_with_padding(A,it)
        szs = collect(size(A))
        strds = collect(strides(A))
        strdsR = collect(strides(ret))

        for dim in 1:N
            n = szs[dim]
            szs[dim] = 1

            M, b = prefiltering_system(eltype(A), n+2pad, it)

            @nloops N i d->1:szs[d] begin
                cc = @ntuple N i
                strt = 1 + sum([(cc[i]-1)*strds[i] for i in 1:length(cc)])
                strtR = 1 + sum([(cc[i]-1)*strdsR[i] for i in 1:length(cc)])
                rng = range(strt, strds[dim], n)
                rngR = range(strtR, strdsR[dim], size(ret,dim))

                b[1+pad:end-pad] += A[rng]
                ret[rngR] = M \ b
                b[1+pad:end-pad], A[rng]
            end
            szs[dim] = n
        end
        ret
    end
end

end # module
