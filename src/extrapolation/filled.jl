nindexes(N::Int) = N == 1 ? "1 index" : "$N indexes"


type FilledInterpolation{T,N,ITP<:AbstractInterpolation,IT,GT,FT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
    fillvalue::FT
end
@doc """
`FilledInterpolation(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.

By comparison with `extrapolate`, this version lets you control the `fillvalue`'s type directly.  It's important for the `fillvalue` to be of the same type as returned by `itp[x1,x2,...]` for in-bounds regions for the index types you are using; otherwise, indexing will be type-unstable (and slow).
""" ->
function FilledInterpolation{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue)
    FilledInterpolation{T,N,typeof(itp),IT,GT,typeof(fillvalue)}(itp, fillvalue)
end

@doc """
`extrapolate(itp, fillvalue)` creates an extrapolation object that returns the `fillvalue` any time the indexes in `itp[x1,x2,...]` are out-of-bounds.
""" ->
extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, fillvalue) = FilledInterpolation(itp, convert(eltype(itp), fillvalue))

@generated function getindex{T,N}(fitp::FilledInterpolation{T,N}, args::Number...)
    n = length(args)
    n == N || return error("Must index $(N)-dimensional interpolation objects with $(nindexes(N))")
    meta = Expr(:meta, :inline)
    quote
        $meta
        # Check to see if we're in the extrapolation region, i.e.,
        # out-of-bounds in an index
        @nexprs $N d->((args[d] < 1 || args[d] > size(fitp.itp, d)) && return fitp.fillvalue)
        # In the interpolation region
        return getindex(fitp.itp,args...)
    end
end

getindex{T}(fitp::FilledInterpolation{T,1}, x::Number, y::Int) = y == 1 ? fitp[x] : throw(BoundsError())
