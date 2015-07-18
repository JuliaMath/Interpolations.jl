immutable Linear <: Degree{1} end

function define_indices_d(::Type{BSpline{Linear}}, d, pad)
    symix, symixp, symfx, symx = symbol("ix_",d), symbol("ixp_",d), symbol("fx_",d), symbol("x_",d)
    quote
        $symix = clamp(floor(Int, real($symx)), 1, size(itp, $d)-1)
        $symixp = $symix + 1
        $symfx = $symx - $symix
    end
end

function coefficients(::Type{BSpline{Linear}}, N, d)
    sym, symp, symfx = symbol(string("c_",d)), symbol(string("cp_",d)), symbol(string("fx_",d))
    quote
        $sym = 1 - $symfx
        $symp = $symfx
    end
end

function gradient_coefficients(::Type{BSpline{Linear}}, d)
    sym, symp, symfx = symbol(string("c_",d)), symbol(string("cp_",d)), symbol(string("fx_",d))
    quote
        $sym = -1
        $symp = 1
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen{IT<:DimSpec{BSpline}}(::Type{BSpline{Linear}}, ::Type{IT}, N::Integer, offsets...)
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
