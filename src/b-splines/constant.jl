struct Constant <: Degree{0} end

"""
Constant b-splines are *nearest-neighbor* interpolations, and effectively
return `A[round(Int,x)]` when interpolating
"""
Constant

function base_rem(::Constant, bounds, x)
    xm = roundbounds(x, bounds)
    δx = x - xm
    fast_trunc(Int, xm), δx
end

expand_index(::Constant, xi::Number, ax::AbstractUnitRange, δx) = (xi,)

value_weights(::Constant, δx) = (oneunit(δx),)
gradient_weights(::Constant, δx) = (zero(δx),)
hessian_weights(::Constant, δx) = (zero(δx),)

padded_axis(ax::AbstractUnitRange, ::BSpline{Constant}) = ax
