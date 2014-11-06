
abstract ExtrapolationBehavior
type ExtrapError <: ExtrapolationBehavior end

function extrap_gen(::Type{ExtrapError}, N)
    quote 
        @nexprs $N d->(1 <= x_d <= size(itp,d) || throw(BoundsError()))
    end
end
