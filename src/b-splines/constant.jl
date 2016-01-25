immutable Constant <: Degree{0} end

"""
Constant b-splines are *nearest-neighbor* interpolations, and effectively
return `A[round(Int,x)]` when interpolating
"""
Constant

"""
`define_indices_d` for a constant b-spline calculates `ix_d = round(Int,x_d)`
"""
function define_indices_d(::Type{BSpline{Constant}}, d, pad)
    symix, symx = symbol("ix_",d), symbol("x_",d)
    :($symix = clamp(round(Int, $symx), 1, size(itp, $d)))
end

"""
`coefficients` for a constant b-spline simply sets `c_d = 1` for compatibility
with the general b-spline framework
"""
function coefficients(::Type{BSpline{Constant}}, N, d)
    sym, symx = symbol("c_",d), symbol("x_",d)
    :($sym = 1)
end

"""
`gradient_coefficients` for a constant b-spline simply sets `c_d = 0` for
compatibility with the general b-spline framework
"""
function gradient_coefficients(::Type{BSpline{Constant}}, d)
    sym, symx = symbol("c_",d), symbol("x_",d)
    :($sym = 0)
end

"""
`hessian_coefficients` for a constant b-spline simply sets `c_d = 0` for
compatibility with the general b-spline framework
"""
function hessian_coefficients(::Type{BSpline{Constant}}, d)
    sym = symbol("c_",d)
    :($sym = 0)
end

function index_gen{IT<:DimSpec{BSpline}}(::Type{BSpline{Constant}}, ::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        d = length(offsets)+1
        sym = symbol("c_", d)
        return :($sym * $(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
