struct Linear{BC<:Union{Throw{OnGrid},Periodic{OnCell}}} <: DegreeBC{1}
    bc::BC
    function Linear{BC}(bc::BC=BC()) where BC<:Union{Throw{OnGrid},Periodic{OnCell}}
        new{BC}(bc)
    end
end

Linear() = Linear(Throw(OnGrid()))
Linear(::Periodic{Nothing}) = Linear(Periodic(OnCell()))
Linear(bc::BC) where BC<:Union{Throw{OnGrid},Periodic{OnCell}} = Linear{BC}(bc)

function Base.show(io::IO, deg::Linear{Throw{OnGrid}})
    print(io, nameof(typeof(deg)), '(', ')')
end

"""
    Linear()

Indicate that the corresponding axis should use linear interpolation.

# Extended help

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

function positions(deg::Linear, ax::AbstractUnitRange{<:Integer}, x)
    x_value = just_dual_value.(x)
    f = floor(x_value)
    # When x == last(ax) we want to use the x-1, x pair
    f = ifelse(x_value == last(ax), f - oneunit(f), f)
    fi = fast_trunc(Int, f)

    expand_index(deg, fi, ax), x - f # for this δ, we want x, not x_value
end
expand_index(::Linear{Throw{OnGrid}}, fi::Number, ax::AbstractUnitRange) = fi
expand_index(::Linear{Periodic{OnCell}}, fi::Number, ax::AbstractUnitRange) =
    (modrange(fi, ax), modrange(fi+1, ax))

value_weights(::Linear, δx) = (1-δx, δx)
gradient_weights(::Linear, δx) = (-oneunit(δx), oneunit(δx))
hessian_weights(::Linear, δx) = (zero(δx), zero(δx))

padded_axis(ax::AbstractUnitRange, ::BSpline{<:Linear}) = ax
