function extrap_prep{N,d}(::Type{Flat}, ::Val{N}, ::Val{d})
    xs_d = Symbol("xs_", d)
    :($xs_d = clamp($xs_d, lbound(etp, $d), ubound(etp, $d)))
end

function extrap_prep{N,d}(::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:lo})
    xs_d = Symbol("xs_", d)
    :($xs_d = max($xs_d, lbound(etp, $d)))
end

function extrap_prep{N,d}(::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:hi})
    xs_d = Symbol("xs_", d)
    :($xs_d = min($xs_d, ubound(etp, $d)))
end

function extrap_prep{N,d}(::Val{:gradient}, ::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:lo})
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]

    quote
        if $xs_d < lbound(etp, $d)
            $xs_d = lbound(etp, $d)
            gradient!(g, etp.itp, $(coords...))
            g[$d] = 0
            return g
        end
    end
end

function extrap_prep{N,d}(::Val{:gradient}, ::Type{Flat}, ::Val{N}, ::Val{d}, ::Val{:hi})
    coords = [Symbol("xs_", k) for k in 1:N]
    xs_d = coords[d]

    quote
        if $xs_d > ubound(etp, $d)
            $xs_d = ubound(etp, $d)
            gradient!(g, etp.itp, $(coords...))
            g[$d] = 0
            return g
        end
    end
end
