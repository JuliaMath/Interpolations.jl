module Interpolations

using Base.Cartesian

import Base: size, eltype, getindex

export BoundaryCondition,
    BCnone,

    ExtrapolationBehavior,
    ExtrapError,
    ExtrapNaN,

    InterpolationDegree,
    Linear,

    AbstractInterpolation,
    Interpolation

abstract BoundaryCondition
immutable BCnone <: BoundaryCondition end

abstract ExtrapolationBehavior
immutable ExtrapError <: ExtrapolationBehavior end
immutable ExtrapNaN <: ExtrapolationBehavior end
immutable ExtrapConstant <: ExtrapolationBehavior end

abstract InterpolationDegree
immutable Linear <: InterpolationDegree end

abstract AbstractInterpolation{T,N,D,BC<:BoundaryCondition,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}

type Interpolation{T,N,ID<:InterpolationDegree,BC<:BoundaryCondition,EB<:ExtrapolationBehavior} <: AbstractInterpolation{T,N,ID,BC,EB}
	coefs::Array{T,N}
end
Interpolation{T,N,ID<:InterpolationDegree,BC<:BoundaryCondition,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Type{ID}, ::Type{BC}, ::Type{EB}) = Interpolation{T,N,ID,BC,EB}(A)

include("linear.jl")

size(itp::Interpolation, d::Integer) = size(itp.coefs, d)
size(itp::Interpolation) = size(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

promote_type_grid(T, x...) = promote_type(T, typeof(x)...)

# This creates getindex methods for all supported combinations
for ID in (Linear,)
    # for BC in subtypes(BoundaryCondition)
    for BC in (BCnone,)
        for EB in (ExtrapError,ExtrapNaN,ExtrapConstant)
            eval(ngenerate(:N, :(promote_type_grid(T, x...)), :(getindex{T,N}(itp::Interpolation{T,N,$ID,$BC,$EB}, x::NTuple{N,Real}...)),
                      N->body_gen(ID, BC, EB, N)))
        end
    end
end


end # module
