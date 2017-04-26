function extrap_prep{N,d}(::Type{Linear}, ::Val{N}, ::Val{d}, ::Val{:lo})
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]
    quote
        if $xs_d < lbound(etp.itp, $d, inds_etp[$d])
            $xs_d = lbound(etp.itp, $d, inds_etp[$d])
            return etp[$(coords...)] + gradient(etp, $(coords...))[$d] * (xs[$d] - $xs_d)
        end
    end
end
function extrap_prep{N,d}(::Type{Linear}, ::Val{N}, ::Val{d}, ::Val{:hi})
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]
    quote
        if $xs_d > ubound(etp, $d, inds_etp[$d])
            $xs_d = ubound(etp, $d, inds_etp[$d])
            return etp[$(coords...)] + gradient(etp, $(coords...))[$d] * (xs[$d] - $xs_d)
        end
    end
end

extrap_prep{N,d,l}(::Val{:gradient}, ::Type{Linear}, n::Val{N}, dim::Val{d}, lohi::Val{l}) =
    extrap_prep(Flat, n, dim, lohi)

extrap_prep{N,d}(::Val{:gradient}, ::Type{Linear}, n::Val{N}, dim::Val{d}) =
    extrap_prep(Flat, n, dim)
