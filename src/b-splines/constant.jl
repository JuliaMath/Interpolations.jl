export Nearest, Previous, Next

abstract type ConstantInterpType end
struct Nearest <: ConstantInterpType end
struct Previous <: ConstantInterpType end
struct Next <: ConstantInterpType end

struct Constant{T<:ConstantInterpType} <: Degree{0}
    Constant() = new{Nearest}()
    Constant{Previous}() = new{Previous}()
    Constant{Next}() = new{Next}()
end

"""
Constant b-splines are *nearest-neighbor* interpolations, and effectively
return `A[round(Int,x)]` when interpolating.
"""
Constant

function positions(c::Constant{Previous}, ax, x)  # discontinuity occurs at integer locations
    xm = floorbounds(x, ax)
    δx = x - xm
    fast_trunc(Int, xm), δx
end
function positions(c::Constant{Next}, ax, x)  # discontinuity occurs at integer locations
    xm = ceilbounds(x, ax)
    δx = x - xm
    fast_trunc(Int, xm), δx
end
function positions(c::Constant{Nearest}, ax, x)  # discontinuity occurs at half-integer locations
    xm = roundbounds(x, ax)
    δx = x - xm
    fast_trunc(Int, xm), δx
end

value_weights(::Constant, δx) = (1,)
gradient_weights(::Constant, δx) = (0,)
hessian_weights(::Constant, δx) = (0,)

padded_axis(ax::AbstractUnitRange, ::BSpline{<:Constant}) = ax
