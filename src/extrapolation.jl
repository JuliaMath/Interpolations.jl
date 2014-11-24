
abstract ExtrapolationBehavior
type ExtrapError <: ExtrapolationBehavior end

function extrap_gen(::OnGrid, ::ExtrapError, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || throw(BoundsError()))
    end
end
extrap_gen(::OnCell, e::ExtrapError, N) = extrap_gen(OnGrid(), e, N)

type ExtrapNaN <: ExtrapolationBehavior end

function extrap_gen(::OnGrid, ::ExtrapNaN, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || return convert(T, NaN))
    end
end
extrap_gen(::OnCell, e::ExtrapNaN, N) = extrap_gen(OnGrid(), e, N)

type ExtrapConstant <: ExtrapolationBehavior end
function extrap_gen(::OnGrid, ::ExtrapConstant, N)
    quote
        @nexprs $N d->(x_d = clamp(x_d, 1, size(itp,d)))
    end
end
extrap_gen(::OnCell, e::ExtrapConstant, N) = extrap_gen(OnGrid(), e, N)
