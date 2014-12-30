module Interpolations

using Base.Cartesian
using Compat

import Base: size, eltype, getindex

export 
    Interpolation,
    Constant,
    Linear,
    Quadratic,
    ExtrapError,
    ExtrapNaN,
    ExtrapConstant,
    OnCell,
    OnGrid,
    ExtendInner,
    Flat

abstract Degree{N}

abstract GridRepresentation
type OnGrid <: GridRepresentation end
type OnCell <: GridRepresentation end

abstract BoundaryCondition
type None <: BoundaryCondition end
type ExtendInner <: BoundaryCondition end
type Flat <: BoundaryCondition end

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

size(itp::Interpolation, d::Integer) = size(itp.coefs, d)
size(itp::Interpolation) = size(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")

gridrepresentation{D,BC,GR<:GridRepresentation}(::InterpolationType{D,BC,GR}) = GR()
degree{D<:Degree,BC,GR}(::InterpolationType{D,BC,GR}) = D()

include("constant.jl")
include("linear.jl")
include("quadratic.jl")


# This creates getindex methods for all supported combinations
for IT in (Constant{OnCell},Linear{OnGrid},Quadratic{ExtendInner,OnCell},Quadratic{Flat,OnCell})
    for EB in (ExtrapError,ExtrapNaN,ExtrapConstant)

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
for IT in (Quadratic{ExtendInner,OnCell},Quadratic{Flat,OnCell})
    @ngenerate N promote_type_grid(T, x...) function prefilter{T,N}(A::Array{T,N},it::IT)
        ret = similar(A)
        szs = collect(size(A))
        strds = collect(strides(A))

        for dim in 1:N
            n = szs[dim]
            szs[dim] = 1

            M = prefiltering_system_matrix(eltype(A), n, it)

            @nloops N i d->1:szs[d] begin
                cc = @ntuple N i
                strt = 1 + sum([(cc[i]-1)*strds[i] for i in 1:length(cc)])
                rng = range(strt, strds[dim], n)
                ret[rng] = M \ vec(A[rng])
            end            
            szs[dim] = n
        end
        ret
    end
end

end # module
