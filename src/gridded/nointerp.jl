function define_indices_d(::Type{Gridded{NoInterp}}, d, pad)
    symix, symx = symbol("ix_",d), symbol("x_",d)
    :($symix = convert(Int, $symx))
end

function coefficients(::Type{Gridded{NoInterp}}, N, d)
    :()
end

function index_gen{IT<:DimSpec{Gridded}}(::Type{Gridded{NoInterp}}, ::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        return :($(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
