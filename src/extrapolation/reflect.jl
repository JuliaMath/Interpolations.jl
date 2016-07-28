"""
`extrap_prep_dim(::Type{Reflect}, d)`

First, translate x into the domain over [lbound, 2(ubound-lbound)) i.e. into twice the size of the domain.
Next, if x is now in the upper part of this ''double-domain´´, reflect over the middle to obtain a new value x' for which f(x') == f(x), but where x' is inside the domain
"""
function extrap_prep_dim(::Type{Reflect}, d)
    xs_d = Symbol("xs_", d)
    quote
        start = lbound(etp.itp, $d)
        width = ubound(etp.itp, $d) - start

        $xs_d = mod($xs_d - start, 2width) + start
        $xs_d > start + width && (xs_d = start + width - $xs_d)
    end
end

extrap_prep{N,d}(::Type{Reflect}, ::Val{N}, ::Val{d}) = extrap_prep_dim(Reflect, d)
function extrap_prep{N,d}(::Type{Reflect}, ::Val{N}, ::Val{d}, ::Val{:lo})
    xs_d = Symbol("xs_", d)
    :($xs_d < lbound(etp.itp, $d) && $(extrap_prep_dim(Reflect, d)))
end
function extrap_prep{N,d}(::Type{Reflect}, ::Val{N}, ::Val{d}, ::Val{:hi})
    xs_d = Symbol("xs_", d)
    :($xs_d > ubound(etp.itp, $d) && $(extrap_prep_dim(Reflect, d)))
end
