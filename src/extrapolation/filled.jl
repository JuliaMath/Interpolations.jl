nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"

mutable struct FilledExtrapolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end

function FilledExtrapolation(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT}
    Te = promote_type(T,typeof(fillvalue))
    FilledExtrapolation{Te,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

Base.parent(A::FilledExtrapolation) = A.itp

"""
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.
"""
extrapolate(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT} = FilledExtrapolation(itp, fillvalue)

@inline function getindex(fitp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, args::Vararg{Number,M}) where {T,N,ITP,IT,GT,FT,M}
    inds, trailing = Base.IteratorsMD.split(args, Val{N})
    @boundscheck all(x->x==1, trailing) || Base.throw_boundserror(fitp, args)
    Tret = typeof(prod(inds) * zero(T))
    checkbounds(Bool, fitp, inds...) && return convert(Tret, fitp.itp[inds...])
    convert(Tret, fitp.fillvalue)
end

@inline Base.checkbounds(::Type{Bool}, A::FilledExtrapolation, I...) = _checkbounds(A, 1, indices(A), I)
@inline _checkbounds(A, d::Int, IA::TT1, I::TT2) where {TT1,TT2} =
    (I[1] >= lbound(A, d, IA[1])) & (I[1] <= ubound(A, d, IA[1])) & _checkbounds(A, d+1, Base.tail(IA), Base.tail(I))
_checkbounds(A, d::Int, ::Tuple{}, ::Tuple{}) = true

getindex(fitp::FilledExtrapolation{T,1}, x::Number, y::Int) where {T} = y == 1 ? fitp[x] : throw(BoundsError())

lbound(etp::FilledExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::FilledExtrapolation, d) = ubound(etp.itp, d)
lbound(etp::FilledExtrapolation, d, inds) = lbound(etp.itp, d, inds)
ubound(etp::FilledExtrapolation, d, inds) = ubound(etp.itp, d, inds)
