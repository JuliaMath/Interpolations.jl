struct Constant <: Degree{0}
    mode
    Constant(mode=:nearest) = new(mode)
end

"""
Constant b-splines are *nearest-neighbor* interpolations, and effectively
return `A[round(Int,x)]` when interpolating.
"""
Constant

function positions(c::Constant, ax, x)  # discontinuity occurs at half-integer locations
    if c.mode == :previous
        xm = floorbounds(x, ax)
    elseif c.mode == :next
        xm = ceilbounds(x, ax)
    else
        xm = roundbounds(x, ax)
    end
    δx = x - xm
    fast_trunc(Int, xm), δx
end

value_weights(::Constant, δx) = (1,)
gradient_weights(::Constant, δx) = (0,)
hessian_weights(::Constant, δx) = (0,)

padded_axis(ax::AbstractUnitRange, ::BSpline{Constant}) = ax
