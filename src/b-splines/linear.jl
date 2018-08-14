struct Linear <: Degree{1} end

"""
Assuming uniform knots with spacing 1, the `i`th peice of linear b-spline
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

function base_rem(::Linear, bounds, x)
    xf = floorbounds(x, bounds)
    xf -= ifelse(xf >= floor(bounds[2]), oneunit(xf), zero(xf))
    δx = x - xf
    fast_trunc(Int, xf), δx
end

expand_index(::Linear, xi::Number, ax::AbstractUnitRange, δx) = (xi, xi+(δx>0))

value_weights(::Linear, δx) = (1-δx, δx)
gradient_weights(::Linear, δx) = (-oneunit(δx), oneunit(δx))
hessian_weights(::Linear, δx) = (zero(δx), zero(δx))

padded_axis(ax::AbstractUnitRange, ::BSpline{Linear}) = ax
