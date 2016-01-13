immutable Cubic{BC<:Flag} <: Degree{3} end
Cubic{BC<:Flag}(::BC) = Cubic{BC}()

"""
Assuming uniform knots with spacing 1, the `i`th piece of cubic spline
implemented here is defined as follows.

    y_i(x) = ϕm p(x-i) + ϕ q(x-i) + ϕp q(1- (x-i)) + ϕpp p(1 - (x-i))

where

    p(☆) = 1/6 * (1-☆)^3
    q(☆) = 2/3 - ☆^2 + 1/2 ☆^3

and the values `ϕX` for `X ⋹ {m, _, p, pp}` are the pre-filtered coefficients.

For future reference, this expands out to the following polynomial:

    y_i(x) = 1/6 cm (1+i-x)^3 + c (2/3 - (x-i)^2 + 1/2 (x-i)^3) +
             cp (2/3 - (1+i-x)^2 + 1/2 (1+i-x)^3) + 1/6 cpp (x-i)^3

When we derive boundary conditions we will use derivatives `y_0'(x)` and
`y_0''(x)`
"""
Cubic

function define_indices_d{BC}(::Type{BSpline{Cubic{BC}}}, d, pad)
    symix, symixm, symixp = symbol("ix_",d), symbol("ixm_",d), symbol("ixp_",d)
    symixpp, symx, symfx = symbol("ixpp_",d), symbol("x_",d), symbol("fx_",d)
    quote
        # ensure that all of ix_d, ixm_d, ixp_d, and ixpp_d are in-bounds no
        # matter the value of pad
        $symix = clamp(floor(Int, $symx), $(2-pad), size(itp,$d)+$(pad-2))
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixm = $symix - 1
        $symixp = $symix + 1
        $symixpp = $symixp + 1
    end
end

function define_indices_d(::Type{BSpline{Cubic{Periodic}}}, d, pad)
    symix, symixm, symixp = symbol("ix_",d), symbol("ixm_",d), symbol("ixp_",d)
    symixpp, symx, symfx = symbol("ixpp_",d), symbol("x_",d), symbol("fx_",d)
    quote
        $symix = clamp(floor(Int, $symx), 1, size(itp,$d))
        $symfx = $symx - $symix
        $symixm = mod1($symix - 1, size(itp,$d))
        $symixp = mod1($symix + 1, size(itp,$d))
        $symixpp = mod1($symix + 2, size(itp,$d))
    end
end

padding{BC<:Flag}(::Type{BSpline{Cubic{BC}}}) = Val{1}()
padding(::Type{BSpline{Cubic{Periodic}}}) = Val{0}()

"""
In this function we assume that `fx_d = x-ix_d` and we produce `cX_d` for
`X ⋹ {m, _, p, pp}` such that

    cm_d  = p(fx_d)
    c_d   = q(fx_d)
    cp_d  = q(1-fx_d)
    cpp_d = p(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Cubic`.
"""
function coefficients{C<:Cubic}(::Type{BSpline{C}}, N, d)
    symm, sym =  symbol("cm_",d), symbol("c_",d)
    symp, sympp = symbol("cp_",d) ,symbol("cpp_",d)
    symfx = symbol("fx_",d)
    symfx_cub = symbol("fx_cub_", d)
    sym_1m_fx_cub = symbol("one_m_fx_cub_", d)
    quote
        $symfx_cub = cub($symfx)
        $sym_1m_fx_cub = cub(1-$symfx)
        $symm = SimpleRatio(1,6)*$sym_1m_fx_cub
        $sym  = SimpleRatio(2,3) - sqr($symfx) + SimpleRatio(1,2)*$symfx_cub
        $symp = SimpleRatio(2,3) - sqr(1-$symfx) + SimpleRatio(1,2)*$sym_1m_fx_cub
        $sympp = SimpleRatio(1,6)*$symfx_cub
    end
end

"""
In this function we assume that `fx_d = x-ix_d` and we produce `cX_d` for
`X ⋹ {m, _, p, pp}` such that

    cm_d  = p'(fx_d)
    c_d   = q'(fx_d)
    cp_d  = q'(1-fx_d)
    cpp_d = p'(1-fx_d)

where `p` and `q` are defined in the docstring for `Cubic`.
"""
function gradient_coefficients{C<:Cubic}(::Type{BSpline{C}}, d)
    symm, sym, symp, sympp = symbol("cm_",d), symbol("c_",d), symbol("cp_",d), symbol("cpp_",d)
    symfx = symbol("fx_",d)
    symfx_sqr = symbol("fx_sqr_", d)
    sym_1m_fx_sqr = symbol("one_m_fx_sqr_", d)
    quote
        $symfx_sqr = sqr($symfx)
        $sym_1m_fx_sqr = sqr(1 - $symfx)

        $symm  = -SimpleRatio(1,2) * $sym_1m_fx_sqr
        $sym   =  SimpleRatio(3,2) * $symfx_sqr     - 2 * $symfx
        $symp  = -SimpleRatio(3,2) * $sym_1m_fx_sqr + 2 * (1 - $symfx)
        $sympp =  SimpleRatio(1,2) * $symfx_sqr
    end
end

function index_gen{C<:Cubic,IT<:DimSpec{BSpline}}(::Type{BSpline{C}}, ::Type{IT}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp, sympp =  symbol("cm_",d), symbol("c_",d), symbol("cp_",d), symbol("cpp_",d)
        return :($symm * $(index_gen(IT, N, offsets...,-1)) + $sym * $(index_gen(IT, N, offsets..., 0)) +
                 $symp * $(index_gen(IT, N, offsets..., 1)) + $sympp * $(index_gen(IT, N, offsets..., 2)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

# ------------ #
# Prefiltering #
# ------------ #

function inner_system_diags{T,C<:Cubic}(::Type{T}, n::Int, ::Type{C})
    du = fill(convert(T, SimpleRatio(1, 6)), n-1)
    d = fill(convert(T, SimpleRatio(2, 3)), n)
    dl = copy(du)
    dl, d, du
end

"""
`Flat` `OnGrid` amounts to setting `y_0'(x) = 0` at `x = 0`. Applying this
condition gives:

    -cm + cp = 0
"""
function prefiltering_system{T,TC}(::Type{T}, ::Type{TC}, n::Int,
                                   ::Type{Cubic{Flat}}, ::Type{OnGrid})
    dl, d, du = inner_system_diags(T, n, Cubic{Flat})
    d[1] = d[end] = -one(T)
    du[1] = dl[end] = zero(T)

    # Now Woodbury correction to set `[1, 3], [n, n-2] ==> 1`
    specs = _build_woodbury_specs(T, n, (1, 3, one(T)), (n, n-2, one(T)))

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), specs...), zeros(TC, n)
end

"""
`Flat` `OnCell` amounts to setting `y_0'(x) = 0` at `x = -1/2`. Applying this
condition gives:

    -9 cm + 11 c -3 cp + 1 cpp = 0
"""
function prefiltering_system{T,TC}(::Type{T}, ::Type{TC}, n::Int,
                                   ::Type{Cubic{Flat}}, ::Type{OnCell})
    dl, d, du = inner_system_diags(T,n,Cubic{Flat})
    d[1] = d[end] = -9
    du[1] = dl[end] = 11

    # now need Woodbury correction to set :
    #    - [1, 3] and [n, n-2] ==> -3
    #    - [1, 4] and [n, n-3] ==> 1
    specs = _build_woodbury_specs(T, n,
                                  (1, 3, T(-3)),
                                  (n, n-2, T(-3)),
                                  (1, 4, one(T)),
                                  (n, n-3, one(T))
                                  )

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), specs...), zeros(TC, n)
end

"""
`Line` `OnCell` amounts to setting `y_0''(x) = 0` at `x = -1/2`. Applying this
condition gives:

    3 cm -7 c + 5 cp -1 cpp = 0
"""
function prefiltering_system{T,TC}(::Type{T}, ::Type{TC}, n::Int,
                                   ::Type{Cubic{Line}}, ::Type{OnCell})
    dl,d,du = inner_system_diags(T,n,Cubic{Flat})
    d[1] = d[end] = 3
    du[1] = dl[end] = -7

    # now need Woodbury correction to set :
    #    - [1, 3] and [n, n-2] ==> -3
    #    - [1, 4] and [n, n-3] ==> 1
    specs = _build_woodbury_specs(T, n,
                                  (1, 3, T(5)),
                                  (n, n-2, T(5)),
                                  (1, 4, -one(T)),
                                  (n, n-3, -one(T))
                                  )

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), specs...), zeros(TC, n)
end

"""
`Line` `OnGrid` amounts to setting `y_0''(x) = 0` at `x = 0`. Applying this
condition gives:

    1 cm -2 c + 1 cp = 0 = 0

This is the same system as `Quadratic{Line}`, `OnGrid` so we reuse the
implementation
"""
function prefiltering_system{T,TC,GT<:GridType}(::Type{T}, ::Type{TC}, n::Int,
                                                ::Type{Cubic{Line}}, ::Type{GT})
    prefiltering_system(T, TC, n, Quadratic{Line}, OnGrid)
end

function prefiltering_system{T,TC,GT<:GridType}(::Type{T}, ::Type{TC}, n::Int,
                                                ::Type{Cubic{Periodic}}, ::Type{GT})
    dl, d, du = inner_system_diags(T,n,Cubic{Periodic})

    specs = _build_woodbury_specs(T, n,
                                  (1, n, SimpleRatio(1, 6)),
                                  (n, 1, SimpleRatio(1, 6))
                                  )

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), specs...), zeros(TC, n)
end

"""
The free boundary condition for either `OnGrid` or `OnCell` makes sure the
interpoland has a continuous third derivative at the second-to-outermost cell
boundary: `y_0'''(1) = y_1'''(1)` and `y_{n-1}'''(n) = y_n'''(n)`. Applying this
condition gives:

    1 cm -3 c + 3 cp -1 cpp = 0

This is the same system as `Quadratic{Free}` so we reuse the implementation
"""
function prefiltering_system{T,TC,GT<:GridType}(::Type{T}, ::Type{TC}, n::Int,
                                                ::Type{Cubic{Free}}, ::Type{GT})
    prefiltering_system(T, TC, n, Quadratic{Free}, GT)
end
