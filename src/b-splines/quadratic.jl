struct Quadratic{BC<:Flag} <: Degree{2} end
Quadratic(::BC) where {BC<:Flag} = Quadratic{BC}()

"""
Assuming uniform knots with spacing 1, the `i`th piece of quadratic spline
implemented here is defined as follows:

    y_i(x) = cm p(x-i) + c q(x) + cp p(1-(x-i))

where

    p(δx) = (δx - 1)^2 / 2
    q(δx) = 3/4 - δx^2

and the values for `cX` for `X ∈ {m,_,p}` are the pre-filtered coefficients.

For future reference, this expands to the following polynomial:

    y_i(x) = cm * 1/2 * (x-i-1)^2 + c * (3/4 - x + i)^2 + cp * 1/2 * (x-i)^2

When we derive boundary conditions we will use derivatives `y_1'(x-1)` and
`y_1''(x-1)`
"""
Quadratic

"""
`define_indices_d` for a quadratic b-spline calculates `ix_d = floor(x_d)` and
`fx_d = x_d - ix_d` (corresponding to `i` `and `δx` in the docstring for
`Quadratic`), as well as auxiliary quantities `ixm_d` and `ixp_d`
"""
function define_indices_d(::Type{BSpline{Quadratic{BC}}}, d, pad) where BC
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symx, symfx = Symbol("x_",d), Symbol("fx_",d)
    quote
        # ensure that all three ix_d, ixm_d, and ixp_d are in-bounds no matter
        # the value of pad
        $symix = clamp(round(Int, $symx), first(inds_itp[$d])+1-$pad, last(inds_itp[$d])+$pad-1)
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixp = $symix + 1
        $symixm = $symix - 1
    end
end
function define_indices_d(::Type{BSpline{Quadratic{Periodic}}}, d, pad)
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symx, symfx = Symbol("x_",d), Symbol("fx_",d)
    quote
        $symix = clamp(round(Int, $symx), first(inds_itp[$d]), last(inds_itp[$d]))
        $symfx = $symx - $symix
        $symixp = modrange($symix + 1, inds_itp[$d])
        $symixm = modrange($symix - 1, inds_itp[$d])
    end
end
function define_indices_d(::Type{BSpline{Quadratic{BC}}}, d, pad) where BC<:Union{InPlace,InPlaceQ}
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symx, symfx = Symbol("x_",d), Symbol("fx_",d)
    pad == 0 || error("Use $BC only with interpolate!")
    quote
        # ensure that all three ix_d, ixm_d, and ixp_d are in-bounds no matter
        # the value of pad
        $symix = clamp(round(Int, $symx), first(inds_itp[$d]), last(inds_itp[$d]))
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixp = min(last(inds_itp[$d]), $symix + 1)
        $symixm = max(first(inds_itp[$d]), $symix - 1)
    end
end

"""
In `coefficients` for a quadratic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p}` such that

    cm_d  = p(fx_d)
    c_d   = q(fx_d)
    cp_d  = p(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Quadratic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function coefficients(::Type{BSpline{Q}}, N, d) where Q<:Quadratic
    symm, sym, symp =  Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d)
    symfx = Symbol("fx_",d)
    quote
        $symm = sqr($symfx - SimpleRatio(1,2))/2
        $sym  = SimpleRatio(3,4) - sqr($symfx)
        $symp = sqr($symfx + SimpleRatio(1,2))/2
    end
end

"""
In `gradient_coefficients` for a quadratic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p}` such that

    cm_d  = p'(fx_d)
    c_d   = q'(fx_d)
    cp_d  = p'(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Quadratic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function gradient_coefficients(::Type{BSpline{Q}}, d) where Q<:Quadratic
    symm, sym, symp =  Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d)
    symfx = Symbol("fx_",d)
    quote
        $symm = $symfx - SimpleRatio(1,2)
        $sym = -2 * $symfx
        $symp = $symfx + SimpleRatio(1,2)
    end
end

"""
In `hessian_coefficients` for a quadratic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p}` such that

    cm_d  = p''(fx_d)
    c_d   = q''(fx_d)
    cp_d  = p''(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Quadratic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function hessian_coefficients(::Type{BSpline{Q}}, d) where Q<:Quadratic
    symm, sym, symp = Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d)
    quote
        $symm = 1
        $sym  = -2
        $symp = 1
    end
end

# This assumes integral values ixm_d, ix_d, and ixp_d,
# coefficients cm_d, c_d, and cp_d, and an array itp.coefs
function index_gen(::Type{BSpline{Q}}, ::Type{IT}, N::Integer, offsets...) where {Q<:Quadratic,IT<:DimSpec{BSpline}}
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp =  Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d)
        return :($symm * $(index_gen(IT, N, offsets...,-1)) + $sym * $(index_gen(IT, N, offsets..., 0)) +
                 $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

function index_gen(::Type{BS}, ::Type{BS}, N::Integer, offsets...) where {BS<:BSpline{<:Quadratic}}
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp =  Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d)
        return :($symm * $(index_gen(BS, BS, N, offsets...,-1)) + $sym * $(index_gen(BS, BS, N, offsets..., 0)) +
                 $symp * $(index_gen(BS, BS, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

padding(::Type{BSpline{Quadratic{BC}}}) where {BC<:Flag} = Val{1}()
padding(::Type{BSpline{Quadratic{Periodic}}}) = Val{0}()

function inner_system_diags(::Type{T}, n::Int, ::Type{Q}) where {T,Q<:Quadratic}
    du = fill(convert(T, SimpleRatio(1,8)), n-1)
    d = fill(convert(T, SimpleRatio(3,4)), n)
    dl = copy(du)
    (dl,d,du)
end

"""
`Quadratic{Flat}` `OnCell` and `Quadratic{Reflect}` `OnCell` amounts to setting
`y_1'(x) = 0` at x=1/2. Applying this condition yields

    -cm + c = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{BC}}, ::Type{OnCell}) where {T,TC,BC<:Union{Flat,Reflect}}
    dl,d,du = inner_system_diags(T,n,Quadratic{BC})
    d[1] = d[end] = -1
    du[1] = dl[end] = 1
    lut!(dl, d, du), zeros(TC, n)
end

function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{InPlace}}, ::Type{OnCell}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,Quadratic{InPlace})
    d[1] = d[end] = convert(T, SimpleRatio(7,8))
    lut!(dl, d, du), zeros(TC, n)
end

# InPlaceQ continues the quadratic at 2 all the way down to 1 (rather than 1.5)
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{InPlaceQ}}, ::Type{OnCell}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,Quadratic{InPlaceQ})
    d[1] = d[end] = SimpleRatio(9,8)
    dl[end] = du[1] = SimpleRatio(-1,4)
    # Woodbury correction to add 1/8 for row 1, col 3 and row n, col n-2
    rowspec = spzeros(T, n, 2)
    colspec = spzeros(T, 2, n)
    valspec = zeros(T, 2, 2)
    valspec[1,1] = valspec[2,2] = SimpleRatio(1,8)
    rowspec[1,1] = rowspec[n,2] = 1
    colspec[1,3] = colspec[2,n-2] = 1
    Woodbury(lut!(dl, d, du), rowspec, valspec, colspec), zeros(TC, n)
end

"""
`Quadratic{Flat}` `OnGrid` and `Quadratic{Reflect}` `OnGrid` amount to setting
`y_1'(x) = 0` at `x=1`. Applying this condition yields

    -cm + cp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{BC}}, ::Type{OnGrid}) where {T,TC,BC<:Union{Flat,Reflect}}
    dl,d,du = inner_system_diags(T,n,Quadratic{BC})
    d[1] = d[end] = -1
    du[1] = dl[end] = 0

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, 3, oneunit(T)),
                                  (n, n-2, oneunit(T))
                                 )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Quadratic{Line}` `OnGrid` and `Quadratic{Line}` `OnCell` amount to setting
`y_1''(x) = 0` at `x=1` and `x=1/2` respectively. Since `y_i''(x)` is independent
of `x` for a quadratic b-spline, these both yield

    1 cm -2 c + 1 cp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{Line}}, ::Type{GT}) where {T,TC,GT<:GridType}
    dl,d,du = inner_system_diags(T,n,Quadratic{Line})
    d[1] = d[end] = 1
    du[1] = dl[end] = -2

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, 3, oneunit(T)),
                                  (n, n-2, oneunit(T)),
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Quadratic{Free}` `OnGrid` and `Quadratic{Free}` `OnCell` amount to requiring
an extra continuous derivative at the second-to-last cell boundary; this means
that `y_1''(3/2) = y_2''(3/2)`, yielding

    1 cm -3 c + 3 cp - cpp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{Free}}, ::Type{GT}) where {T,TC,GT<:GridType}
    dl,d,du = inner_system_diags(T,n,Quadratic{Free})
    d[1] = d[end] = 1
    du[1] = dl[end] = -3

    specs = WoodburyMatrices.sparse_factors(T, n,
                                    (1, 3, 3),
                                    (1, 4, -1),
                                    (n, n-2, 3),
                                    (n, n-3, -1))

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Quadratic{Periodic}` `OnGrid` and `Quadratic{Periodic}` `OnCell` close the system
by looking at the coefficients themselves as periodic, yielding

    c0 = c(N+1)

where `N` is the number of data points.
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, ::Type{Quadratic{Periodic}}, ::Type{GT}) where {T,TC,GT<:GridType}
    dl,d,du = inner_system_diags(T,n,Quadratic{Periodic})

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, n, du[1]),
                                  (n, 1, dl[end])
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end
