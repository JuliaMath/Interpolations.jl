function interpolate(A::AbstractArray, ::NoInterp)
    interpolate(Int, eltype(A), A, NoInterp())
end

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{NoInterp}) = 0

interpdegree(::NoInterp) = NoInterp()

iscomplete(::NoInterp) = true

prefilter(::Type{TWeights}, ::Type{TC}, A::AbstractArray, ::NoInterp) where {TWeights, TC} = A

lbound(ax, ::NoInterp) = first(ax)
ubound(ax, ::NoInterp) = last(ax)

base_rem(::NoInterp, bounds, x::Number) = Int(x), 0

expand_index(::NoInterp, xi::Number, ax::AbstractUnitRange, δx) = (xi,)

value_weights(::NoInterp, δx) = (oneunit(δx),)
gradient_weights(::NoInterp, δx) = (NoInterp(),)
hessian_weights(::NoInterp, δx) = (NoInterp(),)

padded_axis(ax::AbstractUnitRange, ::NoInterp) = ax
