export rrule
"""
    ChainRulesCore.rrule(itp::AbstractInterpolation, x...)

ChainRulesCore.jl `rrule` for integration with automatic differentiation libraries.
Note that it gives the gradient only with respect to the evaluation point `x`,
and not the data inside `itp`.
"""
function ChainRulesCore.rrule(itp::AbstractInterpolation, x...)
    y = itp(x...)
    function interpolate_pullback(Δy)
        nope = ChainRulesCore.@not_implemented "`Interpolations.gradient` does not calculate a gradient with respect to the original data, only the evaluation point"
        (nope, (sum(Δy .* g) for g in Interpolations.gradient(itp, x...))...)
    end
    y, interpolate_pullback
end
