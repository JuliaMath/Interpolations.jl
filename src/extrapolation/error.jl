immutable Throw end
type ErrorExtrapolation{T,N,ITPT,IT,GT} <: AbstractExtrapolation{T,N,ITPT,IT,GT}
    itp::IT
end
ErrorExtrapolation{T,ITPT,IT,GT}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}) =
    ErrorExtrapolation{T,N,ITPT,IT,GT}(itp)
extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Throw}) = 
    ErrorExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,N,ITPT,IT,GT}(etp::Type{ErrorExtrapolation{T,N,ITPT,IT,GT}}, xs...)
    :(@nexprs $N d->(lbound(etp,d) <= xs[d] <= ubound(etp,d) || throw(BoundsError())))
end

lbound(etp::ErrorExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::ErrorExtrapolation, d) = ubound(etp.itp, d)
