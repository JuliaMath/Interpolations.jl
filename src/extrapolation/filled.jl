nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"

type FilledExtrapolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end

function FilledExtrapolation{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue)
    FilledExtrapolation{T,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

"""
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.
"""
extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) = FilledExtrapolation(itp, convert(eltype(itp), fillvalue))

function getindex_impl{T,N,ITP,IT,GT,FT}(fitp::Type{FilledExtrapolation{T,N,ITP,IT,GT,FT}}, args)
   n = length(args)
   n == N || return error("Must index $(N)-dimensional interpolation objects with $(nindexes(N))")

   Tret = FT<:Number ? getindex_return_type(ITP, args) : FT
   meta = Expr(:meta, :inline)
   quote
       $meta
       # Check to see if we're in the extrapolation region, i.e.,
       # out-of-bounds in an index
       @nexprs $N d->((args[d] < lbound(fitp,d) || args[d] > ubound(fitp, d))) && return convert($Tret, fitp.fillvalue)::$Tret
       # In the interpolation region
       return convert($Tret, getindex(fitp.itp,args...))::$Tret
   end
end


@generated function getindex{T,N,ITP,IT,GT,FT}(fitp::FilledExtrapolation{T,N,ITP,IT,GT,FT}, args::Number...)
    getindex_impl(fitp, args)
end

getindex{T}(fitp::FilledExtrapolation{T,1}, x::Number, y::Int) = y == 1 ? fitp[x] : throw(BoundsError())

lbound(etp::FilledExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::FilledExtrapolation, d) = ubound(etp.itp, d)
