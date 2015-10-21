# This assumes integral values ixm_d, ix_d, ixp_d, and ixpp_d (typically ixm_d = ix_d-1, ixp_d = ix_d+1, except at boundaries),
# coefficients cm_d, c_d, cp_d, and cpp_d, and an array itp.coefs
function index_gen(::Type{Cubic}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp, sympp =  symbol("cm_",d), symbol("c_",d), symbol("cp_",d), symbol("cpp_",d)
        return :($symm * $(index_gen(Cubic, N, offsets...,-1)) + $sym * $(index_gen(Cubic, N, offsets..., 0)) +
                 $symp * $(index_gen(Cubic, N, offsets..., 1)) + $sympp * $(index_gen(Cubic, N, offsets..., 2)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
