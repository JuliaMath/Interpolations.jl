export Nearest, Previous, Next

abstract type ConstantInterpType end
struct Nearest <: ConstantInterpType end
struct Previous <: ConstantInterpType end
struct Next <: ConstantInterpType end

struct Constant{T<:ConstantInterpType,BC<:Union{Throw{OnGrid},Periodic{OnCell}}} <: DegreeBC{0}
    bc::BC
end

# Default to Nearest and Throw{OnGrid}
Constant() = Constant{Nearest}()
Constant(bc::BoundaryCondition) = Constant{Nearest}(bc)
Constant{T}() where {T<:ConstantInterpType} = Constant{T,Throw{OnGrid}}(Throw(OnGrid()))
Constant{T}(bc::BC) where {T<:ConstantInterpType,BC<:BoundaryCondition} = Constant{T,BC}(bc)

function Base.show(io::IO, deg::Constant)
    print(io, nameof(typeof(deg)), '{', typeof(deg).parameters[1], '}', '(')
    show(io, deg.bc)
    print(io, ')')
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

function positions(c::Constant{Previous,Periodic{OnCell}}, ax, x)
    # We do not use floorbounds because we do not want to add a half at
    # the lowerbound to round up.
    xm = floor(x)
    δx = x - xm
    modrange(fast_trunc(Int, xm), ax), δx
end
function positions(c::Constant{Next,Periodic{OnCell}}, ax, x)  # discontinuity occurs at integer locations
    xm = ceilbounds(x, ax)
    δx = x - xm
    modrange(fast_trunc(Int, xm), ax), δx
end

value_weights(::Constant, δx) = (1,)
gradient_weights(::Constant, δx) = (0,)
hessian_weights(::Constant, δx) = (0,)

padded_axis(ax::AbstractUnitRange, ::BSpline{<:Constant}) = ax
