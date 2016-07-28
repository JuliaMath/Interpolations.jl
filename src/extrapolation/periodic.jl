"""
`extrap_prep_dim(d, ::Type{Periodic})`

Translate x into the domain [lbound, ubound] my means of `mod()`
"""
function extrap_prep_dim(::Type{Periodic}, d)
    xs_d = Symbol("xs_", d)
    :($xs_d = mod(xs[$d] - lbound(etp.itp, $d), ubound(etp.itp, $d) - lbound(etp.itp, $d)) + lbound(etp.itp, $d))
end

extrap_prep{N,d}(::Type{Periodic}, ::Val{N}, ::Val{d}) = extrap_prep_dim(Periodic, d)
function extrap_prep{N,d}(::Type{Periodic}, ::Val{N}, ::Val{d}, ::Val{:lo})
    xs_d = Symbol("xs_", d)
    quote
        if $xs_d < lbound(etp.itp, $d)
            $(extrap_prep_dim(Periodic, d))
        end
    end
end
function extrap_prep{N,d}(::Type{Periodic}, ::Val{N}, ::Val{d}, ::Val{:hi})
    xs_d = Symbol("xs_", d)
    quote
        if $xs_d > ubound(etp.itp, d)
            $(extrap_prep_dim(Periodic, d))
        end
    end
end
