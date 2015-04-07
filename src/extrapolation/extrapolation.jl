abstract ExtrapolationBehavior

function extrapolate{I<:AbstractInterpolation,E<:ExtrapolationBehavior}(itp::I, ::Type{E})
    error("$E extrapolation undefined for $I")
end
