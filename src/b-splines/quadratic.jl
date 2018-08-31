struct Quadratic{BC<:BoundaryCondition} <: DegreeBC{2}
    bc::BC
end
(deg::Quadratic)(gt::GridType) = Quadratic(deg.bc(gt))


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

function positions(deg::Quadratic, ax, x)
    xm = roundbounds(x, ax)
    δx = x - xm
    expand_index(deg, fast_trunc(Int, xm), ax, δx), δx
end

value_weights(::Quadratic, δx) = (
    sqr(δx - SimpleRatio(1,2))/2,
    SimpleRatio(3,4) - sqr(δx),
    sqr(δx + SimpleRatio(1,2))/2)

gradient_weights(::Quadratic, δx) = (
    δx - SimpleRatio(1,2),
    -2 * δx,
    δx + SimpleRatio(1,2))

hessian_weights(::Quadratic, δx) = (oneunit(δx), -2*oneunit(δx), oneunit(δx))

expand_index(::Quadratic, xi::Number, ax::AbstractUnitRange, δx) = xi-1  # uses WeightedAdjIndex
# Others use WeightedArbIndex
expand_index(::Quadratic{<:Periodic}, xi::Number, ax::AbstractUnitRange, δx) =
    (modrange(xi-1, ax), modrange(xi, ax), modrange(xi+1, ax))
expand_index(::Quadratic{BC}, xi::Number, ax::AbstractUnitRange, δx) where BC<:Union{InPlace,InPlaceQ} =
    (max(xi-1, first(ax)), xi, min(xi+1, last(ax)))

padded_axis(ax::AbstractUnitRange, ::BSpline{<:Quadratic}) = first(ax)-1:last(ax)+1
padded_axis(ax::AbstractUnitRange, ::BSpline{Quadratic{BC}}) where BC<:Union{Periodic,InPlace,InPlaceQ} = ax

# # Due to padding we can extend the bounds
# lbound(ax, ::BSpline{Quadratic{BC}}, ::OnGrid) where BC = first(ax) - 0.5
# ubound(ax, ::BSpline{Quadratic{BC}}, ::OnGrid) where BC = last(ax) + 0.5

function inner_system_diags(::Type{T}, n::Int, ::Quadratic) where {T}
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
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{BC}) where {T,TC,BC<:Union{Flat{OnCell},Reflect{OnCell}}}
    dl,d,du = inner_system_diags(T,n,degree)
    d[1] = d[end] = -1
    du[1] = dl[end] = 1
    lut!(dl, d, du), zeros(TC, n)
end

function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{InPlace{OnCell}}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
    d[1] = d[end] = convert(T, SimpleRatio(7,8))
    lut!(dl, d, du), zeros(TC, n)
end

# InPlaceQ continues the quadratic at 2 all the way down to 1 (rather than 1.5)
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{InPlaceQ{OnCell}}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
    d[1] = d[end] = SimpleRatio(9,8)
    dl[end] = du[1] = SimpleRatio(-1,4)
    # Woodbury correction to add 1/8 for row 1, col 3 and row n, col n-2
    rowspec = spzeros(T, n, 2)
    colspec = spzeros(T, 2, n)
    valspec = zeros(T, 2, 2)
    valspec[1,1] = valspec[2,2] = SimpleRatio(1,8)
    rowspec[1,1] = rowspec[end,2] = 1
    colspec[1,3] = colspec[2,end-2] = 1
    Woodbury(lut!(dl, d, du), rowspec, valspec, colspec), zeros(TC, n)
end

"""
`Quadratic{Flat}` `OnGrid` and `Quadratic{Reflect}` `OnGrid` amount to setting
`y_1'(x) = 0` at `x=1`. Applying this condition yields

    -cm + cp = 0
"""
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{BC}) where {T,TC,BC<:Union{Flat{OnGrid},Reflect{OnGrid}}}
    dl,d,du = inner_system_diags(T,n,degree)
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
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{<:Line}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
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
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{<:Free}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)
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
function prefiltering_system(::Type{T}, ::Type{TC}, n::Int, degree::Quadratic{<:Periodic}) where {T,TC}
    dl,d,du = inner_system_diags(T,n,degree)

    specs = WoodburyMatrices.sparse_factors(T, n,
                                  (1, n, du[1]),
                                  (n, 1, dl[end])
                                  )

    Woodbury(lut!(dl, d, du), specs...), zeros(TC, n)
end
