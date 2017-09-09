function extrap_prep(::Type{Linear}, ::Val{N}, ::Val{d}, ::Val{:lo}) where {N,d}
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]
    quote
        if $xs_d < lbound(etp.itp, $d, inds_etp[$d])
            $xs_d = lbound(etp.itp, $d, inds_etp[$d])
            return etp[$(coords...)] + gradient(etp, $(coords...))[$d] * (xs[$d] - $xs_d)
        end
    end
end
function extrap_prep(::Type{Linear}, ::Val{N}, ::Val{d}, ::Val{:hi}) where {N,d}
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]
    quote
        if $xs_d > ubound(etp, $d, inds_etp[$d])
            $xs_d = ubound(etp, $d, inds_etp[$d])
            return etp[$(coords...)] + gradient(etp, $(coords...))[$d] * (xs[$d] - $xs_d)
        end
    end
end

extrap_prep(::Val{:gradient}, ::Type{Linear}, n::Val{N}, dim::Val{d}, lohi::Val{l}) where {N,d,l} =
    extrap_prep(Flat, n, dim, lohi)

extrap_prep(::Val{:gradient}, ::Type{Linear}, n::Val{N}, dim::Val{d}) where {N,d} =
    extrap_prep(Flat, n, dim)
