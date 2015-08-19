type ConstantExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
ConstantExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    ConstantExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Flat}) =
    ConstantExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,ITP,IT}(exp::Type{ConstantExtrapolation{T,1,ITP,IT,OnGrid}}, x)
    :(x = clamp(x, 1, size(exp,1)))
end
function extrap_prep{T,ITP,IT}(exp::Type{ConstantExtrapolation{T,1,ITP,IT,OnCell}}, x)
    :(x = clamp(x, .5, size(exp,1)+.5))
end
function extrap_prep{T,N,ITP,IT}(exp::Type{ConstantExtrapolation{T,N,ITP,IT,OnGrid}}, xs...)
    :(@nexprs $N d->(xs[d] = clamp(xs[d], 1, size(exp,d))))
end
function extrap_prep{T,N,ITP,IT}(exp::Type{ConstantExtrapolation{T,N,ITP,IT,OnCell}}, xs...)
    :(@nexprs $N d->(xs[d] = clamp(xs[d], .5, size(exp,d)+.5)))
end
