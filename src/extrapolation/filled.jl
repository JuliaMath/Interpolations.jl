mutable struct FilledExtrapolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end

function FilledExtrapolation(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT}
    Te = promote_type(T,typeof(fillvalue))
    FilledExtrapolation{Te,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

Base.parent(A::FilledExtrapolation) = A.itp
etpflag(A::FilledExtrapolation) = A.fillvalue
itpflag(A::FilledExtrapolation) = itpflag(A.itp)

"""
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp(x1,x2,...)` are out-of-bounds.
"""
extrapolate(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT} = FilledExtrapolation(itp, fillvalue)

@inline function (etp::FilledExtrapolation{T,N})(x::Vararg{Number,N}) where {T,N}
    itp = parent(etp)
    Tret = typeof(prod(x) * zero(T))
    if checkbounds(Bool, itp, x...)
        convert(Tret, expand_value(itp, x))
    else
        convert(Tret, etp.fillvalue)
    end
end
@inline function (etp::FilledExtrapolation{T,N})(args::Vararg{Number,M}) where {T,M,N}
    inds, trailing = Base.IteratorsMD.split(args, Val(N))
    @boundscheck all(x->x==1, trailing) || Base.throw_boundserror(etp, args)
    @assert length(inds) == N
    etp(inds...)
end

expand_index_resid_etp(deg, fillvalue, (l, u), x, etp::FilledExtrapolation, xN) =
    (l <= x <= u || Base.throw_boundserror(etp, xN))

# expand_etp_valueE(fv::FT, etp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, x) where {T,N,ITP,IT,GT,FT} = fv
# expand_etp_gradientE(fv::FT, etp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, x) where {T,N,ITP,IT,GT,FT} =
#     zero(SVector{N,FT})
# expand_etp_hessianE(fv::FT, etp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, x) where {T,N,ITP,IT,GT,FT} =
#     zero(Matrix{N,N,FT})
