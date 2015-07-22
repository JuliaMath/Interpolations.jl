export Throw,
       FilledExtrapolation   # for direct control over typeof(fillvalue)

type Extrapolation{T,N,ITPT,IT,GT,ET} <: AbstractInterpolationWrapper{T,N,ITPT,IT,GT}
    itp::ITPT
end
Extrapolation{T,ITPT,IT,GT,ET}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}, ::Type{ET}) =
    Extrapolation{T,N,ITPT,IT,GT,ET}(itp)

"""
`extrapolate(itp, scheme)` adds extrapolation behavior to an interpolation object, according to the provided scheme.

The scheme can take any of these values:

* `Throw` - throws a BoundsError for out-of-bounds indices
* `Flat` - for constant extrapolation, taking the closest in-bounds value
* `Linear` - linear extrapolation (the wrapped interpolation object must support gradient)
* `Reflect` - reflecting extrapolation
* `Periodic` - periodic extrapolation

You can also combine schemes in tuples. For example, the scheme `Tuple{Linear, Flat}` will use linear extrapolation in the first dimension, and constant in the second.

Finally, you can specify different extrapolation behavior in different direction. `Tuple{Tuple{Linear,Flat}, Flat}` will extrapolate linearly in the first dimension if the index is too small, but use constant etrapolation if it is too large, and always use constant extrapolation in the second dimension.
"""
extrapolate{T,N,IT,GT,ET}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{ET}) =
    Extrapolation(T,N,itp,IT,GT,ET)

include("throw.jl")
include("flat.jl")
include("linear.jl")
include("reflect.jl")
include("periodic.jl")

include("extrap_prep.jl")
include("extrap_prep_gradient.jl")

function getindex_impl{T,N,ITPT,IT,GT,ET}(etp::Type{Extrapolation{T,N,ITPT,IT,GT,ET}}, xs...)
    coords = [symbol("xs_",d) for d in 1:N]
    quote
        @nexprs $N d->(xs_d = xs[d])
        $(extrap_prep(ET, Val{N}()))
        etp.itp[$(coords...)]
    end
end

@generated function getindex{T,N,ITPT,IT,GT,ET}(etp::Extrapolation{T,N,ITPT,IT,GT,ET}, xs...)
    getindex_impl(etp, xs...)
end


function gradient!_impl{T,N,ITPT,IT,GT,ET}(g, etp::Type{Extrapolation{T,N,ITPT,IT,GT,ET}}, xs...)
    coords = [symbol("xs_", d) for d in 1:N]
    quote
        @nexprs $N d->(xs_d = xs[d])
        $(extrap_prep(Val{:gradient}(), ET, Val{N}()))
        gradient!(g, etp.itp, $(coords...))
    end
end


@generated function gradient!{T,N,ITPT,IT,GT,ET}(g::AbstractVector, etp::Extrapolation{T,N,ITPT,IT,GT,ET}, xs...)
    gradient!_impl(g, etp, xs...)
end

lbound(etp::Extrapolation, d) = lbound(etp.itp, d)
ubound(etp::Extrapolation, d) = ubound(etp.itp, d)
size(etp::Extrapolation, d) = size(etp.itp, d)

include("filled.jl")
