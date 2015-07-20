type ConstantExtrapolation{T,N,ITP,IT,GT} <: AbstractExtrapolation{T,N,ITP,IT,GT}
    itp::ITP
end
ConstantExtrapolation{T,ITP,IT,GT}(::Type{T}, N, itp::ITP, ::Type{IT}, ::Type{GT}) =
    ConstantExtrapolation{T,N,ITP,IT,GT}(itp)

extrapolate{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}, ::Type{Flat}) =
    ConstantExtrapolation(T,N,itp,IT,GT)

function extrap_prep{T,ITP,IT,GT}(etp::Type{ConstantExtrapolation{T,1,ITP,IT,GT}}, x)
    :(x = clamp(x, lbound(etp,1), ubound(etp,1)))
end
function extrap_prep{T,N,ITP,IT,GT}(etp::Type{ConstantExtrapolation{T,N,ITP,IT,GT}}, xs...)
    :(@nexprs $N d->(xs[d] = clamp(xs[d], lbound(etp,d), ubound(etp,d))))
end

lbound(etp::ConstantExtrapolation, d) = lbound(etp.itp, d)
ubound(etp::ConstantExtrapolation, d) = ubound(etp.itp, d)
