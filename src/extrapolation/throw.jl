function extrap_prep(::Type{Throw}, ::Val{N}, ::Val{d}) where {N,d}
    xsym = Symbol("xs_", d)
    :(lbound(etp, $d, inds_etp[$d]) <= $xsym <= ubound(etp, $d, inds_etp[$d]) || throw(BoundsError()))
end

function extrap_prep(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:lo}) where {N,d}
    xsym = Symbol("xs_", d)
    :(lbound(etp, $d, inds_etp[$d]) <= $xsym || throw(BoundsError()))
end

function extrap_prep(::Type{Throw}, ::Val{N}, ::Val{d}, ::Val{:hi}) where {N,d}
    xsym = Symbol("xs_", d)
    :($xsym <= ubound(etp, $d, inds_etp[$d]) || throw(BoundsError()))
end
