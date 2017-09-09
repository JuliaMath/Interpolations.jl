nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"

mutable struct FilledExtrapolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end

function FilledExtrapolation(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT}
    FilledExtrapolation{T,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

Base.parent(A::FilledExtrapolation) = A.itp

"""
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.
"""
extrapolate(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) where {T,N,IT,GT} = FilledExtrapolation(itp, convert(eltype(itp), fillvalue))

function getindex_impl(fitp::Type{FilledExtrapolation{T,N,ITP,IT,GT,FT}}, args) where {T,N,ITP,IT,GT,FT}
   n = length(args)
   n == N || return error("Must index $(N)-dimensional interpolation objects with $(nindexes(N))")

   Tret = FT<:Number ? getindex_return_type(ITP, args) : FT
   meta = Expr(:meta, :inline)
   quote
       $meta
       # Check to see if we're in the extrapolation region, i.e.,
       # out-of-bounds in an index
       inds_etp = indices(fitp)
       @nexprs $N d->((args[d] < lbound(fitp, d, inds_etp[d]) || args[d] > ubound(fitp, d, inds_etp[d]))) && return convert($Tret, fitp.fillvalue)::$Tret
       # In the interpolation region
       return convert($Tret, getindex(fitp.itp,args...))::$Tret
   end
end


@generated function getindex(fitp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, args::Number...) where {T,N,ITP,IT,GT,FT}
    getindex_impl(fitp, args)
end

getindex(fitp::FilledExtrapolation{T,1}, x::Number, y::Int) where {T} = y == 1 ? fitp[x] : throw(BoundsError())

lbound(etp::FilledExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::FilledExtrapolation, d) = ubound(etp.itp, d)
lbound(etp::FilledExtrapolation, d, inds) = lbound(etp.itp, d, inds)
ubound(etp::FilledExtrapolation, d, inds) = ubound(etp.itp, d, inds)
