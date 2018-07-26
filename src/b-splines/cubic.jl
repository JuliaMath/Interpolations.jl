struct Cubic{BC<:Flag} <: Degree{3} end
Cubic(::BC) where {BC<:Flag} = Cubic{BC}()

"""
Assuming uniform knots with spacing 1, the `i`th piece of cubic spline
implemented here is defined as follows.

    y_i(x) = cm p(x-i) + c q(x-i) + cp q(1- (x-i)) + cpp p(1 - (x-i))

where

    p(δx) = 1/6 * (1-δx)^3
    q(δx) = 2/3 - δx^2 + 1/2 δx^3

and the values `cX` for `X ∈ {m, _, p, pp}` are the pre-filtered coefficients.

For future reference, this expands out to the following polynomial:

    y_i(x) = 1/6 cm (1+i-x)^3 + c (2/3 - (x-i)^2 + 1/2 (x-i)^3) +
             cp (2/3 - (1+i-x)^2 + 1/2 (1+i-x)^3) + 1/6 cpp (x-i)^3

When we derive boundary conditions we will use derivatives `y_0'(x)` and
`y_0''(x)`
"""
Cubic

"""
`define_indices_d` for a cubic b-spline calculates `ix_d = floor(x_d)` and
`fx_d = x_d - ix_d` (corresponding to `i` `and `δx` in the docstring for
`Cubic`), as well as auxiliary quantities `ixm_d`, `ixp_d` and `ixpp_d`
"""
function define_indices_d(::Type{BSpline{Cubic{BC}}}, d, pad) where BC
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symixpp, symx, symfx = Symbol("ixpp_",d), Symbol("x_",d), Symbol("fx_",d)
    quote
        # ensure that all of ix_d, ixm_d, ixp_d, and ixpp_d are in-bounds no
        # matter the value of pad
        $symix = clamp(floor(Int, $symx), first(inds_itp[$d]) + $(1-pad), last(inds_itp[$d]) + $(pad-2))
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixm = $symix - 1
        $symixp = $symix + 1
        $symixpp = $symixp + 1
    end
end

"""
`define_indices_d`  for a cubic, periodic b-spline calculates `ix_d = floor(x_d)`
and `fx_d = x_d - ix_d` (corresponding to `i` and `δx` in the docstring entry
for `Cubic`), as well as auxiliary quantities `ixm_d`, `ixp_d` and `ixpp_d`.

If any `ixX_d` for `x ∈ {m, p, pp}` (note: not `c_d`) should fall outside of
the data interval, they wrap around.
"""
function define_indices_d(::Type{BSpline{Cubic{Periodic}}}, d, pad)
    symix, symixm, symixp = Symbol("ix_",d), Symbol("ixm_",d), Symbol("ixp_",d)
    symixpp, symx, symfx = Symbol("ixpp_",d), Symbol("x_",d), Symbol("fx_",d)
    quote
        tmp = inds_itp[$d]
        $symix = clamp(floor(Int, $symx), first(tmp), last(tmp))
        $symfx = $symx - $symix
        $symixm = modrange($symix - 1, tmp)
        $symixp = modrange($symix + 1, tmp)
        $symixpp = modrange($symix + 2, tmp)
    end
end

padding(::Type{BSpline{Cubic{BC}}}) where {BC<:Flag} = Val{1}()
padding(::Type{BSpline{Cubic{Periodic}}}) = Val{0}()

"""
In `coefficients` for a cubic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p, pp}` such that

    cm_d  = p(fx_d)
    c_d   = q(fx_d)
    cp_d  = q(1-fx_d)
    cpp_d = p(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Cubic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function coefficients(::Type{BSpline{C}}, N, d) where C<:Cubic
    symm, sym =  Symbol("cm_",d), Symbol("c_",d)
    symp, sympp = Symbol("cp_",d) ,Symbol("cpp_",d)
    symfx = Symbol("fx_",d)
    symfx_cub = Symbol("fx_cub_", d)
    sym_1m_fx_cub = Symbol("one_m_fx_cub_", d)
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
In `gradient_coefficients` for a cubic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p, pp}` such that

    cm_d  = p'(fx_d)
    c_d   = q'(fx_d)
    cp_d  = q'(1-fx_d)
    cpp_d = p'(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Cubic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function gradient_coefficients(::Type{BSpline{C}}, d) where C<:Cubic
    symm, sym, symp, sympp = Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d), Symbol("cpp_",d)
    symfx = Symbol("fx_",d)
    symfx_sqr = Symbol("fx_sqr_", d)
    sym_1m_fx_sqr = Symbol("one_m_fx_sqr_", d)
    quote
        $symfx_sqr = sqr($symfx)
        $sym_1m_fx_sqr = sqr(1 - $symfx)

        $symm  = -SimpleRatio(1,2) * $sym_1m_fx_sqr
        $sym   =  SimpleRatio(3,2) * $symfx_sqr     - 2 * $symfx
        $symp  = -SimpleRatio(3,2) * $sym_1m_fx_sqr + 2 * (1 - $symfx)
        $sympp =  SimpleRatio(1,2) * $symfx_sqr
    end
end

"""
In `hessian_coefficients` for a cubic b-spline we assume that `fx_d = x-ix_d`
and we define `cX_d` for `X ⋹ {m, _, p, pp}` such that

    cm_d  = p''(fx_d)
    c_d   = q''(fx_d)
    cp_d  = q''(1-fx_d)
    cpp_d = p''(1-fx_d)

where `p` and `q` are defined in the docstring entry for `Cubic`, and
`fx_d` in the docstring entry for `define_indices_d`.
"""
function hessian_coefficients(::Type{BSpline{C}}, d) where C<:Cubic
    symm, sym, symp, sympp = Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d), Symbol("cpp_",d)
    symfx = Symbol("fx_",d)
    quote
        $symm  = 1 - $symfx
        $sym   = 3 * $symfx - 2
        $symp  = 1 - 3 * $symfx
        $sympp = $symfx
    end
end

function index_gen(::Type{BSpline{C}}, ::Type{IT}, N::Integer, offsets...) where {C<:Cubic,IT<:DimSpec{BSpline}}
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp, sympp =  Symbol("cm_",d), Symbol("c_",d), Symbol("cp_",d), Symbol("cpp_",d)
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

"""
`Cubic`: continuity in function value, first and second derivatives yields

    2/3 1/6
    1/6 2/3 1/6
        1/6 2/3 1/6
           ⋱  ⋱   ⋱
"""
function inner_system_diags(::Type{T}, n::Int, ::Type{C}) where {T,C<:Cubic}
    du = fill(convert(T, SimpleRatio(1, 6)), n-1)
    d = fill(convert(T, SimpleRatio(2, 3)), n)
    dl = copy(du)
    dl, d, du
end

"""
`Cubic{Flat}` `OnGrid` amounts to setting `y_1'(x) = 0` at `x = 1`.
Applying this condition yields

    -cm + cp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Flat}}, ::Type{OnGrid}) where {T,TC}
    dl, d, du = inner_system_diags(T, n, Cubic{Flat})
    d[1] = d[end] = -oneunit(T)
    du[1] = dl[end] = zero(T)

    # Now Woodbury correction to set `[1, 3], [n, n-2] ==> 1`
    specs = WoodburyMatrices.sparse_factors(T, n, (1, 3, oneunit(T)), (n, n-2, oneunit(T)))

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Cubic{Flat}`, `OnCell` amounts to setting `y_1'(x) = 0` at `x = 1/2`.
Applying this condition yields

    -9/8 cm + 11/8 c - 3/8 cp + 1/8 cpp = 0

or, equivalently,

    -9 cm + 11 c -3 cp + 1 cpp = 0

(Note that we use `y_1'(x)` although it is strictly not valid in this domain; if we
were to use `y_0'(x)` we would have to introduce new coefficients, so that would not
close the system. Instead, we extend the outermost polynomial for an extra half-cell.)
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Flat}}, ::Type{OnCell}) where {T,TC}
    dl, d, du = inner_system_diags(T,n,Cubic{Flat})
    d[1] = d[end] = -9
    du[1] = dl[end] = 11

    # now need Woodbury correction to set :
    #    - [1, 3] and [n, n-2] ==> -3
    #    - [1, 4] and [n, n-3] ==> 1
    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, 3, T(-3)),
                                  (n, n-2, T(-3)),
                                  (1, 4, oneunit(T)),
                                  (n, n-3, oneunit(T))
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Cubic{Line}` `OnCell` amounts to setting `y_1''(x) = 0` at `x = 1/2`.
Applying this condition yields

    3 cm -7 c + 5 cp -1 cpp = 0

(Note that we use `y_1'(x)` although it is strictly not valid in this domain; if we
were to use `y_0'(x)` we would have to introduce new coefficients, so that would not
close the system. Instead, we extend the outermost polynomial for an extra half-cell.)
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Line}}, ::Type{OnCell}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,Cubic{Line})
    d[1] = d[end] = 3
    du[1] = dl[end] = -7

    # now need Woodbury correction to set :
    #    - [1, 3] and [n, n-2] ==> -3
    #    - [1, 4] and [n, n-3] ==> 1
    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, 3, T(5)),
                                  (n, n-2, T(5)),
                                  (1, 4, -oneunit(T)),
                                  (n, n-3, -oneunit(T))
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Cubic{Line}` `OnGrid` amounts to setting `y_1''(x) = 0` at `x = 1`. Applying this
condition gives:

    1 cm -2 c + 1 cp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Line}}, ::Type{OnGrid}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,Cubic{Line})
    d[1] = d[end] = 1
    du[1] = dl[end] = -2

    # now need Woodbury correction to set :
    #    - [1, 3] and [n, n-2] ==> 1
    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, 3, oneunit(T)),
                                  (n, n-2, oneunit(T)),
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Cubic{Periodic}` `OnGrid` closes the system by looking at the coefficients themselves
as periodic, yielding

    c0 = c(N+1)

where `N` is the number of data points.
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Periodic}}, ::Type{GT}) where {T,TC,GT<:GridType}
    dl, d, du = inner_system_diags(T,n,Cubic{Periodic})

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, n, du[1]),
                                  (n, 1, dl[end])
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end

"""
`Cubic{Free}` `OnGrid` and `Cubic{Free}` `OnCell` amount to requiring an extra
continuous derivative at the second-to-last cell boundary; this means
`y_1'''(2) = y_2'''(2)`, yielding

    1 cm -3 c + 3 cp -1 cpp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int,
                             ::Type{Cubic{Free}}, ::Type{GT}) where {T,TC,GT<:GridType}
    dl, d, du = inner_system_diags(T,n,Cubic{Periodic})

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, n, du[1]),
                                  (n, 1, dl[end])
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end
