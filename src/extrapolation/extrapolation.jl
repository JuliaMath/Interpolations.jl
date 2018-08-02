mutable struct Extrapolation{T,N,ITPT,IT,GT,ET} <: AbstractExtrapolation{T,N,ITPT,IT,GT}
    itp::ITPT
end

Base.parent(A::Extrapolation) = A.itp

# DimSpec{Flag} is not enough for extrapolation dispatch, since we allow nested tuples
# However, no tuples should be nested deeper than this; the first level is for different
# schemes in different dimensions, and the second level is for different schemes in
# different directions.
const ExtrapDimSpec = Union{Flag,Tuple{Vararg{Union{Flag,NTuple{2,Flag}}}}}

"""
`extrapolate(itp, scheme)` adds extrapolation behavior to an interpolation object, according to the provided scheme.

The scheme can take any of these values:

* `Throw` - throws a BoundsError for out-of-bounds indices
* `Flat` - for constant extrapolation, taking the closest in-bounds value
* `Linear` - linear extrapolation (the wrapped interpolation object must support gradient)
* `Reflect` - reflecting extrapolation (indices must support `mod`)
* `Periodic` - periodic extrapolation (indices must support `mod`)

You can also combine schemes in tuples. For example, the scheme `(Linear(), Flat())` will use linear extrapolation in the first dimension, and constant in the second.

Finally, you can specify different extrapolation behavior in different direction. `((Linear(),Flat()), Flat())` will extrapolate linearly in the first dimension if the index is too small, but use constant etrapolation if it is too large, and always use constant extrapolation in the second dimension.
"""
extrapolate(itp::AbstractInterpolation{T,N,IT,GT}, ::ET) where {T,N,IT,GT,ET<:ExtrapDimSpec} =
    Extrapolation{T,N,typeof(itp),IT,GT,ET}(itp)

count_interp_dims(::Type{<:Extrapolation{T,N,ITPT}}, n) where {T,N,ITPT} = count_interp_dims(ITPT, n)

include("throw.jl")
include("flat.jl")
include("linear.jl")
include("reflect.jl")
include("periodic.jl")

include("extrap_prep.jl")
include("extrap_prep_gradient.jl")

"""
`getindex_impl(::Type{E<:Extrapolation}, xs...)`

Generates an expression to be used
as the function body of the getindex method for the given type of extrapolation
and indices. The heavy lifting is done by the `extrap_prep` function; see
`?extrap_prep` for details.
"""
function getindex_impl(etp::Type{Extrapolation{T,N,ITPT,IT,GT,ET}}, xs...) where {T,N,ITPT,IT,GT,ET}
    coords = [Symbol("xs_",d) for d in 1:N]
    quote
        $(Expr(:meta, :inline))
        @nexprs $N d->(xs_d = xs[d])
        $(extrap_prep(ET, Val{N}()))
        etp.itp[$(coords...)]
    end
end

@generated function getindex(etp::Extrapolation{T,N,ITPT,IT,GT,ET}, xs::Number...) where {T,N,ITPT,IT,GT,ET}
    getindex_impl(etp, xs...)
end

function (etp::Extrapolation{T,N,ITPT,IT,GT,ET})(args...) where {T,N,ITPT,IT,GT,ET}
    # support function calls
    etp[args...]
end

checkbounds(::AbstractExtrapolation,I...) = nothing


function gradient!_impl(g, etp::Type{Extrapolation{T,N,ITPT,IT,GT,ET}}, xs...) where {T,N,ITPT,IT,GT,ET}
    coords = [Symbol("xs_", d) for d in 1:N]
    quote
        $(Expr(:meta, :inline))
        @nexprs $N d->(xs_d = xs[d])
        $(extrap_prep(Val{:gradient}(), ET, Val{N}()))
        gradient!(g, etp.itp, $(coords...))
    end
end


@generated function gradient!(g::AbstractVector, etp::Extrapolation{T,N,ITPT,IT,GT,ET}, xs...) where {T,N,ITPT,IT,GT,ET}
    gradient!_impl(g, etp, xs...)
end

lbound(etp::Extrapolation, d) = lbound(etp.itp, d)
ubound(etp::Extrapolation, d) = ubound(etp.itp, d)
lbound(etp::Extrapolation, d, inds) = lbound(etp.itp, d, inds)
ubound(etp::Extrapolation, d, inds) = ubound(etp.itp, d, inds)
size(etp::Extrapolation, d) = size(etp.itp, d)
@inline axes(etp::AbstractExtrapolation) = axes(etp.itp)
axes(etp::AbstractExtrapolation, d) = axes(etp.itp, d)

include("filled.jl")
