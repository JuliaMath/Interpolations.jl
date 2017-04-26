function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d})
    xsym = Symbol("xs_", d)
    :(lbound(etp, $d, inds_etp[$d]) <= $xsym <= ubound(etp, $d, inds_etp[$d]) || throw(BoundsError()))
end

function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:lo})
    xsym = Symbol("xs_", d)
    :(lbound(etp, $d, inds_etp[$d]) <= $xsym || throw(BoundsError()))
end

function extrap_prep{N,d}(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:hi})
    xsym = Symbol("xs_", d)
    :($xsym <= ubound(etp, $d, inds_etp[$d]) || throw(BoundsError()))
end
