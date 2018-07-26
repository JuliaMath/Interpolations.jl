function define_indices_d(::Type{Gridded{Linear}}, d, pad)
    symix, symixp, symx = Symbol("ix_",d), Symbol("ixp_",d), Symbol("x_",d)
    quote
        $symix = clamp($symix, 1, size(itp, $d)-1)
        $symixp = $symix + 1
    end
end

function coefficients(::Type{Gridded{Linear}}, N, d)
    symix, symixp, symx = Symbol("ix_",d), Symbol("ixp_",d), Symbol("x_",d)
    sym, symp, symfx = Symbol("c_",d), Symbol("cp_",d), Symbol("fx_",d)
    symk, symkix = Symbol("k_",d), Symbol("kix_",d)
    quote
        $symkix = $symk[$symix]
        $symfx = ($symx - $symkix)/($symk[$symixp] - $symkix)
        $sym = 1 - $symfx
        $symp = $symfx
    end
end

function gradient_coefficients(::Type{Gridded{Linear}}, d)
    sym, symp = Symbol("c_",d), Symbol("cp_",d)
    symk, symix = Symbol("k_",d), Symbol("ix_",d)
    symixp = Symbol("ixp_",d)
    quote
        $symp = 1/($symk[$symixp] - $symk[$symix])
        $sym = - $symp
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen(::Type{Gridded{Linear}}, ::Type{IT}, N::Integer, offsets...) where IT<:DimSpec{Gridded}
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
