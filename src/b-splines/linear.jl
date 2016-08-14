immutable Linear <: Degree{1} end

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

"""
`define_indices_d` for a linear b-spline calculates `ix_d = floor(x_d)` and
`fx_d = x_d - ix_d` (corresponding to `i` and `δx` in the docstring for
`Linear`), as well as the auxiliary quantity `ixp_d`
"""
function define_indices_d(::Type{BSpline{Linear}}, d, pad)
    symix, symixp, symfx, symx = Symbol("ix_",d), Symbol("ixp_",d), Symbol("fx_",d), Symbol("x_",d)
    quote
        $symix = clamp(floor(Int, $symx), 1, size(itp, $d)-1)
        $symixp = $symix + 1
        $symfx = $symx - $symix
    end
end

"""
In `coefficients` for a linear b-spline we assume that `fx_d = x-ix-d` and
we define `cX_d` for `X ∈ {_, p}` such that

    c_d = p(fx_d)
    cp_d = p(1-fx_d)

where `p` is defined in the docstring entry for `Linear` and `fx_d` in the
docstring entry for `define_indices_d`.
"""
function coefficients(::Type{BSpline{Linear}}, N, d)
    sym, symp, symfx = Symbol("c_",d), Symbol("cp_",d), Symbol("fx_",d)
    quote
        $sym = 1 - $symfx
        $symp = $symfx
    end
end

"""
In `gradient_coefficients` for a linear b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {_, p}` such that

    c_d  = p'(fx_d)
    cp_d = p'(1-fx_d)

where `p` is defined in the docstring entry for `Linear`, and `fx_d` in the
docstring entry for `define_indices_d`.
"""
function gradient_coefficients(::Type{BSpline{Linear}}, d)
    sym, symp, symfx = Symbol("c_",d), Symbol("cp_",d), Symbol("fx_",d)
    quote
        $sym = -1
        $symp = 1
    end
end

"""
In `hessian_coefficients` for a linear b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {_, p}` such that

    c_d  = p''(fx_d)
    cp_d = p''(1-fx_d)

where `p` is defined in the docstring entry for `Linear`, and `fx_d` in the
docstring entry for `define_indices_d`. (These are both ≡ 0.)
"""
function hessian_coefficients(::Type{BSpline{Linear}}, d)
    sym, symp = Symbol("c_",d), Symbol("cp_",d)
    quote
        $sym = $symp = 0
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen{IT<:DimSpec{BSpline}}(::Type{BSpline{Linear}}, ::Type{IT}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        sym = Symbol("c_", d)
        symp = Symbol("cp_", d)
        return :($sym * $(index_gen(IT, N, offsets..., 0)) + $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
