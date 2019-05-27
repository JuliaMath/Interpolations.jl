struct Cubic{BC<:BoundaryCondition} <: DegreeBC{3}
    bc::BC
end

(deg::Cubic)(gt::GridType) = Cubic(deg.bc(gt))

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

function positions(deg::Cubic, ax, x)
    xf = floorbounds(x, ax)
    xf -= ifelse(xf > last(ax)-1, oneunit(xf), zero(xf))
    δx = x - xf
    expand_index(deg, fast_trunc(Int, xf), ax, δx), δx
end

expand_index(::Cubic{BC}, xi::Number, ax::AbstractUnitRange, δx) where BC = xi-1
expand_index(::Cubic{Periodic{GT}}, xi::Number, ax::AbstractUnitRange, δx) where GT<:GridType =
    (modrange(xi-1, ax), modrange(xi, ax), modrange(xi+1, ax), modrange(xi+2, ax))

function value_weights(::Cubic, δx)
    x3, xcomp3 = cub(δx), cub(1-δx)
    (SimpleRatio(1,6) * xcomp3,
     SimpleRatio(2,3) - sqr(δx) + SimpleRatio(1,2)*x3,
     SimpleRatio(2,3) - sqr(1-δx) + SimpleRatio(1,2)*xcomp3,
     SimpleRatio(1,6) * x3)
end

function gradient_weights(::Cubic, δx)
    x2, xcomp2 = sqr(δx), sqr(1-δx)
    (-SimpleRatio(1,2) * xcomp2,
     -2*δx + SimpleRatio(3,2)*x2,
     +2*(1-δx) - SimpleRatio(3,2)*xcomp2,
     SimpleRatio(1,2) * x2)
end

hessian_weights(::Cubic, δx) = (1-δx, 3*δx-2, 3*(1-δx)-2, δx)


# ------------ #
# Prefiltering #
# ------------ #

padded_axis(ax::AbstractUnitRange, ::BSpline{<:Cubic}) = first(ax)-1:last(ax)+1
padded_axis(ax::AbstractUnitRange, ::BSpline{Cubic{Periodic{GT}}}) where GT<:GridType = ax

# # Due to padding we can extend the bounds
# lbound(ax, ::BSpline{Cubic{BC}}, ::OnGrid) where BC = first(ax) - 0.5
# ubound(ax, ::BSpline{Cubic{BC}}, ::OnGrid) where BC = last(ax) + 0.5

"""
`Cubic`: continuity in function value, first and second derivatives yields

    2/3 1/6
    1/6 2/3 1/6
        1/6 2/3 1/6
           ⋱  ⋱   ⋱
"""
function inner_system_diags(::Type{T}, n::Int, ::Cubic) where {T}
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
                             degree::Cubic{Flat{OnGrid}}) where {T,TC}
    dl, d, du = inner_system_diags(T, n, degree)
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
                             degree::Cubic{Flat{OnCell}}) where {T,TC}
    dl, d, du = inner_system_diags(T,n,degree)
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
                             degree::Cubic{Line{OnCell}}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
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
                             degree::Cubic{Line{OnGrid}}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
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
                             degree::Cubic{<:Periodic}) where {T,TC}
    dl, d, du = inner_system_diags(T,n,degree)

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
                             degree::Cubic{<:Free}) where {T,TC}
    dl, d, du = inner_system_diags(T,n,degree)

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, n, du[1]),
                                  (n, 1, dl[end])
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end
