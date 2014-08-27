module Interpolations

using Base.Cartesian

import Base: size, eltype, getindex

export BoundaryCondition,
    BCnone,

    ExtrapolationBehavior,
    ExtrapError,
    ExtrapNaN,

    InterpolationDegree,
    LinearInterpolation,

    AbstractInterpolation,
    Interpolation

abstract BoundaryCondition
type BCnone <: BoundaryCondition end
abstract ExtrapolationBehavior
type ExtrapError <: ExtrapolationBehavior end
type ExtrapNaN <: ExtrapolationBehavior end
abstract InterpolationDegree
type LinearInterpolation <: InterpolationDegree end

abstract AbstractInterpolation{T,N,D,BC<:BoundaryCondition,EB<:ExtrapolationBehavior} <: AbstractArray{T,N}

type Interpolation{T,N,ID<:InterpolationDegree,BC<:BoundaryCondition,EB<:ExtrapolationBehavior} <: AbstractInterpolation{T,N,ID,BC,EB}
	coefs::Array{T,N}
end
Interpolation{T,N,ID<:InterpolationDegree,BC<:BoundaryCondition,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Type{ID}, ::Type{BC}, ::Type{EB}) = Interpolation{T,N,ID,BC,EB}(A)
# Linear interpolation needs no BC
Interpolation{T,N,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Type{LinearInterpolation}, ::Type{EB}) = Interpolation(A, LinearInterpolation, BCnone, EB)
# typealias InterpLinear{T,N,BC,EB} Interpolation{T,N,1,BC,EB}
# InterpLinear{T,N,BC,EB}(A::Array{T,N}, ::Type{BC}, ::Type{EB}) = Interpolation(A, 1, BC, EB)

size(itp::Interpolation, d::Integer) = size(itp.coefs, d)
size(itp::Interpolation) = size(itp.coefs)
eltype(itp::Interpolation) = eltype(itp.coefs)

function extrap_gen(::Type{ExtrapError}, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || throw(BoundsError()))
    end
end

function extrap_gen(::Type{ExtrapNaN}, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || return(nan(eltype(itp))))
    end
end

function interp_gen(::Type{LinearInterpolation}, N)
    quote
        # ix_d is the index in dimension d of the nearest node *before* the interpolation point
        # ixp_d is the index in dimension d of the nearest node *after* the interpolation point
        # x_d is the d:th coordinate of the interpolation point
        # fx_d is a parameter in [0,1] such that x_d = ix_d + fx_d
        # use @nexprs to generate expressions that set these in all N dimensions
        @nexprs $N d->(ix_d = ifloor(x_d); fx_d = x_d - convert(typeof(x_d), ix_d))
        @nexprs $N d->(ixp_d = ix_d + 1)
        @inbounds ret = $(index_gen(LinearInterpolation, N))
        ret
    end
end

function index_gen(::Type{LinearInterpolation}, N::Integer, offsets...)
    if length(offsets) < N
        sym = symbol("fx_"*string(length(offsets)+1))
        return :((one($sym)-$sym) * $(index_gen(LinearInterpolation, N, offsets..., 0)) + $sym * $(index_gen(LinearInterpolation, N, offsets..., 1)))
    else
        indices = [offsets[d] == 0 ? symbol("ix_"*string(d)) : symbol("ixp_"*string(d)) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

function body_gen{BC,EB}(::Type{LinearInterpolation}, ::Type{BC}, ::Type{EB}, N::Integer)
    extrap_ex = extrap_gen(EB,N)
    interp_ex = interp_gen(LinearInterpolation, N)
    quote
        $extrap_ex
        $interp_ex
    end
end

promote_type_grid(T, x...) = promote_type(T, typeof(x)...)

# This creates getindex methods for all supported combinations
for ID in (LinearInterpolation,)
    # for BC in subtypes(BoundaryCondition)
    for BC in (BCnone,)
        for EB in (ExtrapError,ExtrapNaN)
            eval(ngenerate(:N, :(promote_type_grid(T, x...)), :(getindex{T,N}(itp::Interpolation{T,N,$ID,$BC,$EB}, x::NTuple{N,Real}...)),
                      N->body_gen(ID, BC, EB, N)))
        end
    end
end


end # module
