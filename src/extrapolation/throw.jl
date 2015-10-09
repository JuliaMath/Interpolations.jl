function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d})
    xsym = symbol("xs_", d)
    :(lbound(etp, $d) <= $xsym <= ubound(etp, $d) || throw(BoundsError()))
end

function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:lo})
    xsym = symbol("xs_", d)
    :(lbound(etp, $d) <= $xsym || throw(BoundsError()))
end

function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:hi})
    xsym = symbol("xs_", d)
    :($xsym <= ubound(etp, $d) || throw(BoundsError()))
end
