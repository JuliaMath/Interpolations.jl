struct Constant <: Degree{0} end

"""
Constant b-splines are *nearest-neighbor* interpolations, and effectively
return `A[round(Int,x)]` when interpolating.
"""
Constant

function positions(::Constant, ax, x)  # discontinuity occurs at half-integer locations
    xm = roundbounds(x, ax)
    δx = x - xm
    fast_trunc(Int, xm), δx
end

value_weights(::Constant, δx) = (1,)
gradient_weights(::Constant, δx) = (0,)
hessian_weights(::Constant, δx) = (0,)

padded_axis(ax::AbstractUnitRange, ::BSpline{Constant}) = ax
