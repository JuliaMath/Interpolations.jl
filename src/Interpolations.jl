module Interpolations

using Base.Cartesian

import Base: size, eltype, getindex

export 
    Interpolation, 
    Linear, 
    ExtrapError

abstract Degree{N}

abstract BoundaryCondition
type BCnone <: BoundaryCondition end

abstract GridRepresentation
type OnGrid <: GridRepresentation end

abstract InterpolationType{D<:Degree,BC<:BoundaryCondition,GR<:GridRepresentation}

include("extrapolation.jl")

abstract AbstractInterpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}
type Interpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior} <: AbstractInterpolation{T,N,IT,EB}
    coefs::Array{T,N}
end
Interpolation{T,N,IT<:InterpolationType,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Type{IT}, ::Type{EB}) = Interpolation{T,N,IT,EB}(A)

size(itp::Interpolation, d::Integer) = size(itp.coefs, d)
size(itp::Interpolation) = size(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

include("linear.jl")

promote_type_grid(T, x...) = promote_type(T, typeof(x)...)

# This creates getindex methods for all supported combinations
for IT in (LinearOnGrid,)
    for EB in (ExtrapError,)
        eval(ngenerate(
            :N,
            :(promote_type_grid(T, x...)),
            :(getindex{T,N}(itp::Interpolation{T,N,$IT,$EB}, x::NTuple{N,Real}...)), 
            N->:($(extrap_gen(EB,N)); $(interp_gen(IT, N)))
        ))
    end
end

end # module
