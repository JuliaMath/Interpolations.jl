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

weightedindex_parts(fs, ::NoInterp, ax, x::Number) = Int(x)

# positions(::NoInterp, ax, x) = (Int(x),), 0

# value_weights(::NoInterp, δx) = (oneunit(δx),)
# gradient_weights(::NoInterp, δx) = (NoInterp(),)
# hessian_weights(::NoInterp, δx) = (NoInterp(),)

padded_axis(ax::AbstractUnitRange, ::NoInterp) = ax
