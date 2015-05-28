immutable Throw end
type ErrorExtrapolation{T,N,ITPT,IT,GT} <: AbstractExtrapolation{T,N,ITPT,IT,GT}
    itp::IT
end
ErrorExtrapolation{T,ITPT,IT,GT}(::Type{T}, N, itp::ITPT, ::Type{IT}, ::Type{GT}) =
    ErrorExtrapolation{T,N,ITPT,IT,GT}(itp)
extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Throw}) = 
    ErrorExtrapolation(T,N,itp,IT,GT)


function extrap_prep{T,N,ITPT,IT}(exp::Type{ErrorExtrapolation{T,N,ITPT,IT,OnGrid}}, xs...)
    :(@nexprs $N d->(@show 1 <= xs[d] <= size(exp,d) || throw(BoundsError())))
end
