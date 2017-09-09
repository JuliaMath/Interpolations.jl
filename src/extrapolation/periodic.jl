"""
`extrap_prep_dim(d, ::Type{Periodic})`

Translate x into the domain [lbound, ubound] my means of `mod()`
"""
function extrap_prep_dim(::Type{Periodic}, d)
    xs_d = Symbol("xs_", d)
    :($xs_d = mod(xs[$d] - lbound(etp.itp, $d, inds_etp[$d]), ubound(etp.itp, $d, inds_etp[$d]) - lbound(etp.itp, $d, inds_etp[$d])) + lbound(etp.itp, $d, inds_etp[$d]))
end

extrap_prep(::Type{Periodic}, ::Val{N}, ::Val{d}) where {N,d} = extrap_prep_dim(Periodic, d)
function extrap_prep(::Type{Periodic}, ::Val{N}, ::Val{d}, ::Val{:lo}) where {N,d}
    xs_d = Symbol("xs_", d)
    quote
        if $xs_d < lbound(etp.itp, $d, inds_etp[$d])
            $(extrap_prep_dim(Periodic, d))
        end
    end
end
function extrap_prep(::Type{Periodic}, ::Val{N}, ::Val{d}, ::Val{:hi}) where {N,d}
    xs_d = Symbol("xs_", d)
    quote
        if $xs_d > ubound(etp.itp, d, inds_etp[$d])
            $(extrap_prep_dim(Periodic, d))
        end
    end
end
