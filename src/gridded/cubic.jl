"""
`define_indices_d`  for a cubic, periodic b-spline calculates `ix_d = floor(x_d)`
and `fx_d = x_d - ix_d` (corresponding to `i` and `δx` in the docstring entry
for `Cubic`), as well as auxiliary quantities `ixm_d`, `ixp_d` and `ixpp_d`.

If any `ixX_d` for `x ∈ {m, p, pp}` (note: not `c_d`) should fall outside of
the data interval, they wrap around.
"""
function define_indices_d(::Type{Gridded{Cubic{Line}}}, d, pad)
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symixpp, symx, symfx = Symbol("ixpp_",d), Symbol("x_",d), Symbol("fx_",d)
    quote
        tmp = inds_itp[$d]
        $symix = clamp(floor(Int, $symx), first(tmp), last(tmp))
        $symfx = $symx - $symix
        $symixm = modrange($symix - 1, tmp)
        $symixp = modrange($symix + 1, tmp)
        $symixpp = modrange($symix + 2, tmp)
    end
end

# function coefficients{C<:Cubic}(::Type{Gridded{C}}, N, d)
#     sym_a, sym_b = Symbol("a_", d), Symbol("b_", d)
#     sym_c, sym_d = Symbol("b_", d), Symbol("c_", d)
#     sym_cp, sym_ap = Symbol("cp_", d), Symbol("cp_", d)
#     sym_h, sym_hp = Symbol("h_", d), Symbol("hp_", d)
#
#
#     # symfx = Symbol("fx_",d)
#     # symfx_cub = Symbol("fx_cub_", d)
#     # sym_1m_fx_cub = Symbol("one_m_fx_cub_", d)
#     quote
#         $sym_b = ($sym_ap - $sym_a)
#         $sym_d = ($sym_cp - $sym_c) / (2 * $sym_h)
#     end
# end


function coefficients(::Type{Gridded{Cubic{Line}}}, N, d)
    symix, symx = Symbol("ix_",d), Symbol("x_",d)
    symk, symkix = Symbol("k_",d), Symbol("kix_",d)
    symfx = Symbol("fx_",d)
    symfx_sq = Symbol("fx_sq_", d)
    symfx_cub = Symbol("fx_cub_", d)
    quote
        $symkix = $symk[$symix]
        $symfx = ($symx - $symkix)
        $symfx_sq = $symfx * $symfx
        $symfx_cub = $symfx_sq * $symfx
    end
end

function index_gen{IT<:DimSpec{Gridded}}(::Type{Gridded{Cubic{Line}}}, ::Type{IT}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        sym = Symbol("c_", d)
        symp = Symbol("cp_", d)
            return :($sym * $(index_gen(IT, N, offsets..., 0)) + $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        d = length(offsets)
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        sym_a = Symbol("a_", d); sym_ap = Symbol("ap_", d)
        sym_c = Symbol("c_", d); sym_cp = Symbol("cp_", d)
        sym_b = Symbol("b_", d)
        sym_d = Symbol("d_", d)
        sym_h = Symbol("h_", d)
        sym_k = Symbol("k_", d)
        symix = Symbol("ix_", d); symixp = Symbol("ixp_", d)

        symfx = Symbol("fx_",d)
        symfx_sq = Symbol("fx_sq_", d)
        symfx_cub = Symbol("fx_cub_", d)
        return quote
            $sym_a = itp.coefs[$(indices...)]
            $sym_c = itp.coefs2[$(indices...)]
            $sym_h = $sym_k[$symixp] - $sym_k[$symix]
            $sym_b = ($sym_ap - $sym_a)/$sym_h + $sym_h*(2*$sym_c + $sym_cp)/3
            $sym_d = ($sym_cp - $sym_c) / (3 * $sym_h)
            $sym_a + $sym_b * $symfx + $sym_c * $symfx_sq + $sym_d * $symfx_cub
        end
    end
end
