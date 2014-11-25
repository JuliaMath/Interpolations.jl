
abstract ExtrapolationBehavior
type ExtrapError <: ExtrapolationBehavior end

function extrap_gen(::OnGrid, ::ExtrapError, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || throw(BoundsError()))
    end
end
function extrap_gen(::OnCell, ::ExtrapError, N)
    quote
        @nexprs $N d->(.5 <= x_d <= size(itp,d)+.5 || throw(BoundsError()))
    end
end

type ExtrapNaN <: ExtrapolationBehavior end

function extrap_gen(::OnGrid, ::ExtrapNaN, N)
    quote
        @nexprs $N d->(1 <= x_d <= size(itp,d) || return convert(T, NaN))
    end
end
function extrap_gen(::OnCell, ::ExtrapNaN, N)
    quote
        @nexprs $N d->(.5 <= x_d <: size(itp,d)+.5 || return convert(T, NaN))
    end
end

type ExtrapConstant <: ExtrapolationBehavior end
function extrap_gen(::OnGrid, ::ExtrapConstant, N)
    quote
        @nexprs $N d->(x_d = clamp(x_d, 1, size(itp,d)))
    end
end
function extrap_gen(::OnCell, ::ExtrapConstant, N)
    quote
        @nexprs $N d->(x_d = clamp(x_d, .5, size(itp,d)+.5))
    end
end

