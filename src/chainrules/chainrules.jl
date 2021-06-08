export rrule
"""
    ChainRulesCore.rrule(itp::AbstractInterpolation, x...)
ChainRulesCore.jl `rrule` for integration with automatic differentiation libraries.
"""
function ChainRulesCore.rrule(itp::AbstractInterpolation, x...)
    y = itp(x...)
    function pullback(Δy)
        (ChainRulesCore.NoTangent(), Δy * Interpolations.gradient(itp, x...))
    end
    y, pullback
end