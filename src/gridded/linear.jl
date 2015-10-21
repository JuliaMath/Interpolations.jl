function define_indices_d(::Type{Gridded{Linear}}, d, pad)
    symix, symixp, symx = symbol("ix_",d), symbol("ixp_",d), symbol("x_",d)
    quote
        $symix = clamp($symix, 1, size(itp, $d)-1)
        $symixp = $symix + 1
    end
end

function coefficients(::Type{Gridded{Linear}}, N, d)
    symix, symixp, symx = symbol("ix_",d), symbol("ixp_",d), symbol("x_",d)
    sym, symp, symfx = symbol("c_",d), symbol("cp_",d), symbol("fx_",d)
    symk, symkix = symbol("k_",d), symbol("kix_",d)
    quote
        $symkix = $symk[$symix]
        $symfx = ($symx - $symkix)/($symk[$symixp] - $symkix)
        $sym = 1 - $symfx
        $symp = $symfx
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen{IT<:DimSpec{Gridded}}(::Type{Gridded{Linear}}, ::Type{IT}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        sym = symbol("c_"*string(d))
        symp = symbol("cp_"*string(d))
        return :($sym * $(index_gen(IT, N, offsets..., 0)) + $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
