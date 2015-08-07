function define_indices_d(::Type{Gridded{Constant}}, d, pad)
    symix, symx = symbol("ix_",d), symbol("x_",d)
    symk, symkix = symbol("k_",d), symbol("kix_",d)
    # Choose the closer of k[ix] and k[ix+1]
    quote
        $symix = clamp($symix, 1, size(itp, $d)-1)
        $symkix = $symk[$symix]
        cmp = $symk[$symix+1]
        $symix = abs($symx-$symkix) < abs($symx-cmp) ? $symix : $symix+1
    end
end

function coefficients(::Type{Gridded{Constant}}, N, d)
    sym, symx = symbol(string("c_",d)), symbol(string("x_",d))
    :($sym = 1)
end

function index_gen{IT<:DimSpec{Gridded}}(::Type{Gridded{Constant}}, ::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        d = length(offsets)+1
        sym = symbol("c_"*string(d))
        return :($sym * $(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
