
#=
    Interpolate using cubic Hermite splines. The breakpoints in arrays xbp and ybp are assumed to be sorted.
    Evaluate the function in all points of the array xeval.
    Methods:
        "Linear"                yuck
        "FiniteDifference"      classic cubic interpolation, no tension parameter
                                Finite difference can overshoot for non-monotonic data
        "Cardinal"              cubic cardinal splines, uses tension parameter which must be between [0,1]
                                cubin cardinal splines can overshoot for non-monotonic data
                                (increasing tension decreases overshoot)
        "Akima"                 monotonic - tangents are determined at each given point locally,
                                the curve obtained is close to a manually drawn curve, can overshoot for non-monotonic data
        "FritschCarlson"        monotonic - tangents are first initialized, then adjusted if they are not monotonic
                                can overshoot for non-monotonic data
        "FritschButland"        monotonic - faster algorithm (only requires one pass) but somewhat higher apparent "tension"
        "Steffen"               monotonic - also only one pass, results usually between FritschCarlson and FritschButland
    Sources:
        Akima (1970), "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures", doi:10.1145/321607.321609
        Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021.
        Fritsch & Butland (1984), "A Method for Constructing Local Monotone Piecewise Cubic Interpolants", doi:10.1137/0905021.
        Steffen (1990), "A Simple Method for Monotonic Interpolation in One Dimension", http://adsabs.harvard.edu/abs/1990A%26A...239..443S

    Implementation based on http://bl.ocks.org/niclasmattsson/7bceb05fba6c71c78d507adae3d29417
=#

export
    LinearMonotonicInterpolation,
    FiniteDifferenceMonotonicInterpolation,
    CardinalMonotonicInterpolation,
    AkimaMonotonicInterpolation,
    FritschCarlsonMonotonicInterpolation,
    FritschButlandMonotonicInterpolation,
    SteffenMonotonicInterpolation

"""
    MonotonicInterpolationType

Abstract class for all types of monotonic interpolation.
"""
abstract type MonotonicInterpolationType <: InterpolationType end

"""
    LinearMonotonicInterpolation

Simple linear interpolation.
"""
struct LinearMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    FiniteDifferenceMonotonicInterpolation

Classic cubic interpolation, no tension parameter.
Finite difference can overshoot for non-monotonic data.
"""
struct FiniteDifferenceMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    CardinalMonotonicInterpolation(tension)

Cubic cardinal splines, uses `tension` parameter which must be between [0,1]
Cubin cardinal splines can overshoot for non-monotonic data
(increasing tension reduces overshoot).
"""
struct CardinalMonotonicInterpolation{TTension<:Number} <: MonotonicInterpolationType
    tension :: TTension # must be in [0, 1]
end

"""
    AkimaMonotonicInterpolation

Monotonic interpolation based on [Akima (1970)
"A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures",
doi:10.1145/321607.321609](@cite Akima1970).

Tangents are determined at each given point locally,
results are close to manual drawn curves
"""
struct AkimaMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    FritschCarlsonMonotonicInterpolation

Monotonic interpolation based on [Fritsch & Carlson (1980),
"Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021](@cite Fritsch1980).

Tangents are first initialized, then adjusted if they are not monotonic
can overshoot for non-monotonic data
"""
struct FritschCarlsonMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    FritschButlandMonotonicInterpolation

Monotonic interpolation based on [Fritsch & Butland (1984),
"A Method for Constructing Local Monotone Piecewise Cubic Interpolants",
doi:10.1137/0905021](@cite Fritsch1984).

Faster than FritschCarlsonMonotonicInterpolation (only requires one pass)
but somewhat higher apparent "tension".
"""
struct FritschButlandMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    SteffenMonotonicInterpolation

Monotonic interpolation based on [Steffen (1990),
"A Simple Method for Monotonic Interpolation in One Dimension",
http://adsabs.harvard.edu/abs/1990A%26A...239..443S](@cite Steffen1990)

Only one pass, results usually between FritschCarlson and FritschButland.
"""
struct SteffenMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    MonotonicInterpolation

Monotonic interpolation up to third order represented by type, knots and
coefficients.
"""
struct MonotonicInterpolation{T, TCoeffs1, TCoeffs2, TCoeffs3, TEl, TInterpolationType<:MonotonicInterpolationType,
    TKnots<:AbstractVector{<:Number}, TACoeff <: AbstractArray{TEl,1}} <: AbstractInterpolation{T,1,DimSpec{TInterpolationType}}

    it::TInterpolationType
    knots::TKnots
    A::TACoeff # constant parts of piecewise polynomials
    m::Vector{TCoeffs1} # coefficients of linear parts of piecewise polynomials
    c::Vector{TCoeffs2} # coefficients of quadratic parts of piecewise polynomials
    d::Vector{TCoeffs3} # coefficients of cubic parts of piecewise polynomials
end

function Base.:(==)(o1::MonotonicInterpolation, o2::MonotonicInterpolation)
    o1.it == o2.it &&
    o1.knots == o2.knots &&
    o1.A == o2.A &&
    o1.m == o2.m &&
    o1.c == o2.c &&
    o1.d == o2.d
end

size(A::MonotonicInterpolation) = size(A.knots)
axes(A::MonotonicInterpolation) = axes(A.knots)

itpflag(A::MonotonicInterpolation) = A.it
coefficients(A::MonotonicInterpolation) = A.A

function MonotonicInterpolation(
        ::Type{TWeights}, it::TInterpolationType, knots::TKnots, A::AbstractArray{TEl,1},
        m::Vector{TCoeffs1}, c::Vector{TCoeffs2}, d::Vector{TCoeffs3}
    ) where {TWeights, TCoeffs1, TCoeffs2, TCoeffs3, TEl, TInterpolationType<:MonotonicInterpolationType, TKnots<:AbstractVector{<:Number}}

    isconcretetype(TInterpolationType) || error("The spline type must be a leaf type (was $TInterpolationType)")
    isconcretetype(tcoef(A)) || @warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    check_monotonic(knots, A)

    cZero = zero(TWeights)
    if isempty(A)
        T = Base.promote_op(*, typeof(cZero), eltype(A))
    else
        T = typeof(cZero * first(A))
    end

    MonotonicInterpolation{T, TCoeffs1, TCoeffs2, TCoeffs3, TEl, TInterpolationType, TKnots, typeof(A)}(it, knots, A, m, c, d)
end

function interpolate(
        ::Type{TWeights}, ::Type{TCoeffs1}, ::Type{TCoeffs2}, ::Type{TCoeffs3},
        knots::TKnots, A::AbstractArray{TEl,1}, it::TInterpolationType
    ) where {TWeights,TCoeffs1,TCoeffs2,TCoeffs3,TEl,TKnots<:AbstractVector{<:Number},TInterpolationType<:MonotonicInterpolationType}

    check_monotonic(knots, A)

    # first we need to determine tangents (m)
    n = length(knots)
    m, Δ = calcTangents(TCoeffs1, knots, A, it)
    c = Vector{TCoeffs2}(undef, n-1)
    d = Vector{TCoeffs3}(undef, n-1)
    for k ∈ eachindex(c)
        if TInterpolationType == LinearMonotonicInterpolation
            c[k] = zero(TCoeffs2)
            d[k] = zero(TCoeffs3)
        else
            xdiff = knots[k+1] - knots[k]
            c[k] = (3*Δ[k] - 2*m[k] - m[k+1]) / xdiff
            d[k] = (m[k] + m[k+1] - 2*Δ[k]) / (xdiff * xdiff)
        end
    end

    MonotonicInterpolation(TWeights, it, knots, A, m, c, d)
end

function interpolate(
        knots::AbstractVector{<:Number}, A::AbstractArray{TEl,1},
        it::TInterpolationType
    ) where {TEl,TInterpolationType<:MonotonicInterpolationType}

    interpolate(tweight(A),typeof(oneunit(eltype(A)) / oneunit(eltype(knots))),
        typeof(oneunit(eltype(A)) / oneunit(eltype(knots))^2),
        typeof(oneunit(eltype(A)) / oneunit(eltype(knots))^3),knots,A,it)
end

function (itp::MonotonicInterpolation)(x::Number)
    x_value = just_dual_value.(x)
    @boundscheck (checkbounds(Bool, itp, x_value) || Base.throw_boundserror(itp, (x_value,)))
    k = searchsortedfirst(itp.knots, x_value)
    if k > 1
        k -= 1
    end
    xdiff = x - itp.knots[k]
    return itp.A[k] + itp.m[k]*xdiff + itp.c[k]*xdiff*xdiff + itp.d[k]*xdiff*xdiff*xdiff
end

function gradient(itp::MonotonicInterpolation, x::Number)
    return SVector(gradient1(itp, x))
end

function gradient1(itp::MonotonicInterpolation, x::Number)
    @boundscheck (checkbounds(Bool, itp, x) || Base.throw_boundserror(itp, (x,)))
    k = searchsortedfirst(itp.knots, x)
    if k > 1
        k -= 1
    end
    xdiff = x - itp.knots[k]
    return itp.m[k] + 2*itp.c[k]*xdiff + 3*itp.d[k]*xdiff*xdiff
end

function hessian(itp::MonotonicInterpolation, x::Number)
    return SVector(hessian1(itp, x))
end

function hessian1(itp::MonotonicInterpolation, x::Number)
    @boundscheck (checkbounds(Bool, itp, x) || Base.throw_boundserror(itp, (x,)))
    k = searchsortedfirst(itp.knots, x)
    if k > 1
        k -= 1
    end
    xdiff = x - itp.knots[k]
    return 2*itp.c[k] + 6*itp.d[k]*xdiff
end

@inline function check_monotonic(knots, A)
    axes(knots) == axes(A) || throw(DimensionMismatch("knot vector must have the same axes as the corresponding array"))
    issorted(knots) || error("knot-vector must be sorted in increasing order")
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::LinearMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        m[k] = Δ[k]
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::FiniteDifferenceMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        if k == 1   # left endpoint
            m[k] = Δ[k]
        else
            m[k] = (Δ[k-1] + Δ[k]) / 2
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::CardinalMonotonicInterpolation{TTension}) where {TTension, TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        if k == 1   # left endpoint
            m[k] = Δ[k]
        else
            m[k] = (oneunit(TTension) - method.tension) * (y[k+1] - y[k-1]) / (x[k+1] - x[k-1])
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::AkimaMonotonicInterpolation) where {TCoeffs, TEl}

    # based on Akima (1970),
    # "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures",
    # doi:10.1145/321607.321609.

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
    end
    Γ = [3Δ[1] - 2Δ[2]; 2Δ[1] - Δ[2]; Δ; 2Δ[n-1] - Δ[n-2]; 3Δ[n-1] - 2Δ[n-2]]
    for k ∈ eachindex(m)
        δ = abs(Γ[k+3] - Γ[k+2]) + abs(Γ[k+1] - Γ[k])
        if δ > zero(δ)
            α = abs(Γ[k+1] - Γ[k]) / δ
            m[k] = (1-α) * Γ[k+1] + α * Γ[k+2]
        else
            m[k] = 0.5 * Γ[k+1] + 0.5 * Γ[k+2]
        end
    end
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::FritschCarlsonMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    buff = Vector{typeof(first(Δ)^2)}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        if k == 1   # left endpoint
            m[k] = Δ[k]
        else
            # If any consecutive secant lines change sign (i.e. curve changes direction), initialize the tangent to zero.
            # This is needed to make the interpolation monotonic. Otherwise set tangent to the average of the secants.
            if Δ[k-1] * Δ[k] <= zero(Δ[k]^2)
                m[k] = zero(TCoeffs)
            else
                m[k] =  (Δ[k-1] + Δ[k]) / 2.0
            end
        end
    end
    m[n] = Δ[n-1]
    #=
    Fritsch & Carlson derived necessary and sufficient conditions for monotonicity in their 1980 paper. Splines will be
    monotonic if all tangents are in a certain region of the alpha-beta plane, with alpha and beta as defined below.
    A robust choice is to put alpha & beta within a circle around origo with radius 3. The FritschCarlson algorithm
    makes simple initial estimates of tangents and then does another pass over data points to move any outlier tangents
    into the monotonic region. FritschButland & Steffen algorithms make more elaborate first estimates of tangents that
    are guaranteed to lie in the monotonic region, so no second pass is necessary.
    =#

    # Second pass of FritschCarlson: adjust any non-monotonic tangents.
    for k ∈ eachindex(Δ)
        if Δ[k] == zero(TCoeffs)
            m[k] = zero(TCoeffs)
            m[k+1] = zero(TCoeffs)
            continue
        end
        α = m[k] / Δ[k]
        β = m[k+1] / Δ[k]
        τ = 3.0 * oneunit(α) / sqrt(α^2 + β^2)
        if τ < 1.0   # if we're outside the circle with radius 3 then move onto the circle
            m[k] = τ * α * Δ[k]
            m[k+1] = τ * β * Δ[k]
        end
    end
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y :: AbstractVector{TEl}, method :: FritschButlandMonotonicInterpolation) where {TCoeffs, TEl}

    # based on Fritsch & Butland (1984),
    # "A Method for Constructing Local Monotone Piecewise Cubic Interpolants",
    # doi:10.1137/0905021.

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        if k == 1   # left endpoint
            m[k] = Δ[k]
        elseif Δ[k-1] * Δ[k] <= zero(Δ[k]^2)
            m[k] = zero(TCoeffs)
        else
            α = (1.0 + (x[k+1] - x[k]) / (x[k+1] - x[k-1])) / 3.0
            m[k] = Δ[k-1] * Δ[k] / (α*Δ[k] + (1.0 - α)*Δ[k-1])
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::SteffenMonotonicInterpolation) where {TCoeffs, TEl}

    # Steffen (1990),
    # "A Simple Method for Monotonic Interpolation in One Dimension",
    # http://adsabs.harvard.edu/abs/1990A%26A...239..443S

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k ∈ eachindex(Δ)
        Δ[k] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        if k == 1   # left endpoint
            m[k] = Δ[k]
        else
            p = ((x[k+1] - x[k]) * Δ[k-1] + (x[k] - x[k-1]) * Δ[k]) / (x[k+1] - x[k-1])
            m[k] = (sign(Δ[k-1]) + sign(Δ[k])) *
                min(abs(Δ[k-1]), abs(Δ[k]), 0.5*abs(p))
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{<:MonotonicInterpolationType}) = 1

lbounds(itp::MonotonicInterpolation) = (first(itp.knots),)
ubounds(itp::MonotonicInterpolation) = (last(itp.knots),)
