immutable Constant <: Degree{0} end

function define_indices_d(::Type{BSpline{Constant}}, d, pad)
    symix, symx = symbol("ix_",d), symbol("x_",d)
    :($symix = clamp(round(Int, real($symx)), 1, size(itp, $d)))
end

function coefficients(::Type{BSpline{Constant}}, N, d)
    sym, symx = symbol("c_",d), symbol("x_",d)
    :($sym = 1)
end

function gradient_coefficients(::Type{BSpline{Constant}}, d)
    sym, symx = symbol("c_",d), symbol("x_",d)
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
