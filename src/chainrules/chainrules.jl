export rrule
"""
    ChainRulesCore.rrule(itp::AbstractInterpolation, x...)
ChainRulesCore.jl `rrule` for integration with automatic differentiation libraries.
"""
function ChainRulesCore.rrule(itp::AbstractInterpolation, x...)
    y = itp(x...)
    function pullback(Δy)
        (NO_FIELDS, Δy * Interpolations.gradient(itp, x...))
    end
    y, pullback
end