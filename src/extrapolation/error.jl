immutable ErrorExtrapolation <: ExtrapolationBehavior
type ErrorExtrapolation{T,N,IT,GT} <: AbstracExtrapolation{T,N,IT,GT}
    itp::IT
end

stagedfunction getindex(exp::ErrorExtrapolation, )

