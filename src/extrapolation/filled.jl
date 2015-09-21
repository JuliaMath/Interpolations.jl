nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"


type FilledExtrapolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end
"""
`FilledExtrapolation(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.

By comparison with `extrapolate`, this version lets you control the `fillvalue`'s type directly.  It's important for the `fillvalue` to be of the same type as returned by `itp[x1,x2,...]` for in-bounds regions for the index types you are using; otherwise, indexing will be type-unstable (and slow).
"""
function FilledExtrapolation{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue)
    FilledExtrapolation{T,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

"""
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.
"""
extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) = FilledExtrapolation(itp, convert(eltype(itp), fillvalue))

@generated function getindex{T,N}(fitp::FilledExtrapolation{T,N}, args::Number...)
    n = length(args)
    n == N || return error("Must index $(N)-dimensional interpolation objects with $(nindexes(N))")
    meta = Expr(:meta, :inline)
    quote
        $meta
        # Check to see if we're in the extrapolation region, i.e.,
        # out-of-bounds in an index
        @nexprs $N d->((args[d] < lbound(fitp,d) || args[d] > ubound(fitp, d)) && return fitp.fillvalue)
        # In the interpolation region
        return getindex(fitp.itp,args...)
    end
end

getindex{T}(fitp::FilledExtrapolation{T,1}, x::Number, y::Int) = y == 1 ? fitp[x] : throw(BoundsError())

lbound(etp::FilledExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::FilledExtrapolation, d) = ubound(etp.itp, d)
