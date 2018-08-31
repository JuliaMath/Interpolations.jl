struct Linear <: Degree{1} end  # boundary conditions not supported

"""
Assuming uniform knots with spacing 1, the `i`th piece of linear b-spline
implemented here is defined as follows.

    y_i(x) = c p(x) + cp p(1-x)

where

    p(δx) = x

and the values `cX` for `X ∈ {_, p}` are the coefficients.

Linear b-splines are naturally interpolating, and require no prefiltering;
there is therefore no need for boundary conditions to be provided.

Also, although the implementation is slightly different in order to re-use
the framework built for general b-splines, the resulting interpolant is just
a piecewise linear function connecting each pair of neighboring data points.
"""
Linear

function positions(::Linear, ax::AbstractUnitRange{<:Integer}, x)
    f = floor(x)
    # When x == last(ax) we want to use the x-1, x pair
    f = ifelse(x == last(ax), f - oneunit(f), f)
    fi = fast_trunc(Int, f)
    return fi, x-f
end

value_weights(::Linear, δx) = (1-δx, δx)
gradient_weights(::Linear, δx) = (-oneunit(δx), oneunit(δx))
hessian_weights(::Linear, δx) = (zero(δx), zero(δx))

padded_axis(ax::AbstractUnitRange, ::BSpline{Linear}) = ax
