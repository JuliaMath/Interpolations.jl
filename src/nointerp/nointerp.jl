function interpolate(A::AbstractArray, ::NoInterp, gt::GT) where {GT<:DimSpec{GridType}}
    interpolate(Int, eltype(A), A, NoInterp(), gt)
end

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{NoInterp}) = 0

interpdegree(::NoInterp) = NoInterp()

prefilter(::Type{TWeights}, ::Type{TC}, A::AbstractArray, ::NoInterp, ::GridType) where {TWeights, TC} = A

lbound(ax, ::NoInterp, ::GridType) = first(ax)
ubound(ax, ::NoInterp, ::GridType) = last(ax)

base_rem(::NoInterp, bounds, x::Number) = Int(x), 0

expand_index(::NoInterp, xi::Number, ax::AbstractUnitRange, δx) = (xi,)

value_weights(::NoInterp, δx) = (oneunit(δx),)
gradient_weights(::NoInterp, δx) = (zero(δx),)
hessian_weights(::NoInterp, δx) = (zero(δx),)

padded_axis(ax::AbstractUnitRange, ::NoInterp) = ax
