function extrap_prep(::Type{Flat}, ::Val{N}, ::Val{d}) where {N,d}
    xs_d = Symbol("xs_", d)
    :($xs_d = clamp($xs_d, lbound(etp, $d, inds_etp[$d]), ubound(etp, $d, inds_etp[$d])))
end

function extrap_prep(::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:lo}) where {N,d}
    xs_d = Symbol("xs_", d)
    :($xs_d = max($xs_d, lbound(etp, $d, inds_etp[$d])))
end

function extrap_prep(::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:hi}) where {N,d}
    xs_d = Symbol("xs_", d)
    :($xs_d = min($xs_d, ubound(etp, $d, inds_etp[$d])))
end

function extrap_prep(::Val{:gradient}, ::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:lo}) where {N,d}
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]

    quote
        if $xs_d < lbound(etp, $d, inds_etp[$d])
            $xs_d = lbound(etp, $d, inds_etp[$d])
            gradient!(g, etp.itp, $(coords...))
            g[$d] = 0
            return g
        end
    end
end

function extrap_prep(::Val{:gradient}, ::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:hi}) where {N,d}
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]

    quote
        if $xs_d > ubound(etp, $d, inds_etp[$d])
            $xs_d = ubound(etp, $d, inds_etp[$d])
            gradient!(g, etp.itp, $(coords...))
            g[$d] = 0
            return g
        end
    end
end
