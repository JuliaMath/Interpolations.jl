function base_rem(::Constant, knotv, ki, x)
    l, u = knotv[ki], knotv[ki+1]

    xm = roundbounds(x, bounds)
    δx = x - xm
    fast_trunc(Int, xm), δx
end

function define_indices_d(::Type{Gridded{Constant}}, d, pad)
    symix, symx = Symbol("ix_",d), Symbol("x_",d)
    symk, symkix = Symbol("k_",d), Symbol("kix_",d)
    # Choose the closer of k[ix] and k[ix+1]
    quote
        $symix = clamp($symix, 1, size(itp, $d)-1)
        $symkix = $symk[$symix]
        cmp = $symk[$symix+1]
        $symix = abs($symx-$symkix) < abs($symx-cmp) ? $symix : $symix+1
    end
end

function coefficients(::Type{Gridded{Constant}}, N, d)
    sym, symx = Symbol("c_",d), Symbol("x_",d)
    :($sym = 1)
end

function gradient_coefficients(::Type{Gridded{Constant}}, N, d)
    sym, symx = Symbol("c_",d), Symbol("x_",d)
    :($sym = 0)
end

function index_gen(::Type{Gridded{Constant}}, ::Type{IT}, N::Integer, offsets...) where IT<:DimSpec{Gridded}
    if (length(offsets) < N)
        d = length(offsets)+1
        sym = Symbol("c_", d)
        return :($sym * $(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
