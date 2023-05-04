import Adapt: adapt_structure
using Adapt: adapt

function adapt_structure(to, itp::BSplineInterpolation{T,N}) where {T,N}
    coefs′ = adapt(to, itp.coefs)
    T′ = update_eltype(T, coefs′, itp.coefs)
    BSplineInterpolation{T′,N}(coefs′, itp.parentaxes, itp.it)
end

function update_eltype(T, coefs′, coefs)
    ET = eltype(coefs′)
    ET === eltype(coefs) && return T
    WT = tweight(coefs′)
    T′ = Base.promote_op(*, WT, ET)
    (isconcretetype(T′) || isempty(coefs)) && return T′
    return typeof(zero(WT) * convert(ET, first(coefs)))
end

function adapt_structure(to, itp::LanczosInterpolation{T,N}) where {T,N}
    coefs′ = adapt(to, itp.coefs)
    parentaxes′ = adapt(to, itp.parentaxes)
    LanczosInterpolation{eltype(coefs′),N}(coefs′, parentaxes′, itp.it)
end

function adapt_structure(to, itp::GriddedInterpolation{T,N}) where {T,N}
    coefs′ = adapt(to, itp.coefs)
    knots′ = adapt(to, itp.knots)
    T′ = update_eltype(T, coefs′, itp.coefs)
    GriddedInterpolation{T′,N,typeof(coefs′),itptype(itp),typeof(knots′)}(knots′, coefs′, itp.it)
end

function adapt_structure(to, itp::ScaledInterpolation{T,N,<:Any,IT,RT}) where {T,N,IT,RT<:NTuple{N,AbstractRange}}
    ranges = itp.ranges
    itp′ = adapt(to, itp.itp)
    ScaledInterpolation{eltype(itp′),N,typeof(itp′),IT,RT}(itp′, ranges)
end

function adapt_structure(to, itp::Extrapolation{T,N}) where {T,N}
    et = itp.et
    itp′ = adapt(to, itp.itp)
    Extrapolation{eltype(itp′),N,typeof(itp′),itptype(itp),typeof(et)}(itp′, et)
end

function adapt_structure(to, itp::FilledExtrapolation{T,N}) where {T,N}
    fillvalue = itp.fillvalue
    itp′ = adapt(to, itp.itp)
    FilledExtrapolation{eltype(itp′),N,typeof(itp′),itptype(itp),typeof(fillvalue)}(itp′, fillvalue)
end

import Base.Broadcast: broadcasted, BroadcastStyle
using Base.Broadcast: broadcastable, combine_styles, AbstractArrayStyle
function broadcasted(itp::AbstractInterpolation, args...)
    args′ = map(broadcastable, args)
    # we overload BroadcastStyle here (try our best to do broadcast on GPU)
    style = combine_styles(Ref(itp), args′...)
    broadcasted(style, itp, args′...)
end

"""
    Interpolations.root_storage_type(::Type{<:AbstractInterpolation}) -> Type{<:AbstractArray}

This function returns the type of the root coefficients array of an `AbstractInterpolation`.
Some array wrappers, like `OffsetArray`, should be skipped.
"""
root_storage_type(::Type{T}) where {T<:AbstractInterpolation} = Array{eltype(T),ndims(T)} # fallback to `Array` by default.
root_storage_type(::Type{T}) where {T<:Extrapolation} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:FilledExtrapolation} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:ScaledInterpolation} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:BSplineInterpolation} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:LanczosInterpolation} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:GriddedInterpolation} = root_storage_type(fieldtype(T, 2))
root_storage_type(::Type{T}) where {T<:OffsetArray} = root_storage_type(fieldtype(T, 1))
root_storage_type(::Type{T}) where {T<:AbstractArray} = T

BroadcastStyle(::Type{<:Ref{T}}) where {T<:AbstractInterpolation} = _to_scalar_style(BroadcastStyle(T))
BroadcastStyle(::Type{T}) where {T<:AbstractInterpolation} = BroadcastStyle(root_storage_type(T))

_to_scalar_style(::S) where {S<:AbstractArrayStyle} = S(Val(0))
_to_scalar_style(S::AbstractArrayStyle{Any}) = S
_to_scalar_style(S) = S
