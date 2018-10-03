
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
        "FritschCarlson"        monotonic - tangents are first initialized, then adjusted if they are not monotonic
                                can overshoot for non-monotonic data
        "FritschButland"        monotonic - faster algorithm (only requires one pass) but somewhat higher apparent "tension"
        "Steffen"               monotonic - also only one pass, results usually between FritschCarlson and FritschButland
    Sources:
        Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021.
        Fritsch & Butland (1984), "A Method for Constructing Local Monotone Piecewise Cubic Interpolants", doi:10.1137/0905021.
        Steffen (1990), "A Simple Method for Monotonic Interpolation in One Dimension", http://adsabs.harvard.edu/abs/1990A%26A...239..443S

    Implementation based on http://bl.ocks.org/niclasmattsson/7bceb05fba6c71c78d507adae3d29417
=#

export
    LinearMonotonicInterpolation,
    FiniteDifferenceMonotonicInterpolation,
    CardinalMonotonicInterpolation,
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
    FritschCarlsonMonotonicInterpolation

Monotonic interpolation based on Fritsch & Carlson (1980),
"Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021.

Tangents are first initialized, then adjusted if they are not monotonic
can overshoot for non-monotonic data
"""
struct FritschCarlsonMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    FritschButlandMonotonicInterpolation

Monotonic interpolation based on  Fritsch & Butland (1984),
"A Method for Constructing Local Monotone Piecewise Cubic Interpolants",
doi:10.1137/0905021.

Faster than FritschCarlsonMonotonicInterpolation (only requires one pass)
but somewhat higher apparent "tension".
"""
struct FritschButlandMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    SteffenMonotonicInterpolation

Monotonic interpolation based on Steffen (1990),
"A Simple Method for Monotonic Interpolation in One Dimension",
http://adsabs.harvard.edu/abs/1990A%26A...239..443S

Only one pass, results usually between FritschCarlson and FritschButland.
"""
struct SteffenMonotonicInterpolation <: MonotonicInterpolationType
end

"""
    MonotonicInterpolation

Monotonic interpolation up to third order represented by type, knots and
coefficients.
"""
struct MonotonicInterpolation{T, TCoeffs, TEl, TInterpolationType<:MonotonicInterpolationType,
    TKnots<:AbstractVector{<:Number}, TACoeff <: AbstractArray{TEl,1}} <: AbstractInterpolation{T,1,DimSpec{TInterpolationType}}

    it::TInterpolationType
    knots::TKnots
    A::TACoeff # constant parts of piecewise polynomials
    m::Vector{TCoeffs} # coefficients of linear parts of piecewise polynomials
    c::Vector{TCoeffs} # coefficients of quadratic parts of piecewise polynomials
    d::Vector{TCoeffs} # coefficients of cubic parts of piecewise polynomials
end


size(A::MonotonicInterpolation) = size(A.knots)
axes(A::MonotonicInterpolation) = axes(A.knots)

itpflag(A::MonotonicInterpolation) = A.it

function MonotonicInterpolation(::Type{TWeights}, it::TInterpolationType, knots::TKnots, A::AbstractArray{TEl,1},
    m::Vector{TCoeffs}, c::Vector{TCoeffs}, d::Vector{TCoeffs}) where {TWeights, TCoeffs, TEl, TInterpolationType<:MonotonicInterpolationType, TKnots<:AbstractVector{<:Number}}

    isconcretetype(TInterpolationType) || error("The spline type must be a leaf type (was $TInterpolationType)")
    isconcretetype(TCoeffs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    check_monotonic(knots, A)

    cZero = zero(TWeights)
    if isempty(A)
        T = Base.promote_op(*, typeof(cZero), eltype(A))
    else
        T = typeof(cZero * first(A))
    end

    MonotonicInterpolation{T, TCoeffs, TEl, TInterpolationType, TKnots, typeof(A)}(it, knots, A, m, c, d)
end

function interpolate(::Type{TWeights}, ::Type{TCoeffs}, knots::TKnots,
    A::AbstractArray{TEl,1}, it::TInterpolationType) where {TWeights,TCoeffs,TEl,TKnots<:AbstractVector{<:Number},TInterpolationType<:MonotonicInterpolationType}

    check_monotonic(knots, A)

    # first we need to determine tangents (m)
    n = length(knots)
    m, Δ = calcTangents(TCoeffs, knots, A, it)
    c = Vector{TCoeffs}(undef, n-1)
    d = Vector{TCoeffs}(undef, n-1)
    for k ∈ 1:n-1
        if TInterpolationType == LinearMonotonicInterpolation
            c[k] = d[k] = zero(TCoeffs)
        else
            xdiff = knots[k+1] - knots[k]
            c[k] = (3*Δ[k] - 2*m[k] - m[k+1]) / xdiff
            d[k] = (m[k] + m[k+1] - 2*Δ[k]) / (xdiff * xdiff)
        end
    end

    MonotonicInterpolation(TWeights, it, knots, A, m, c, d)
end

function interpolate(knots::AbstractVector{<:Number}, A::AbstractArray{TEl,1},
    it::TInterpolationType) where {TEl,TInterpolationType<:MonotonicInterpolationType}

    interpolate(tweight(A), tcoef(A), knots, A, it)
end

function (itp::MonotonicInterpolation)(x::Number)
    @boundscheck (checkbounds(Bool, itp, x) || Base.throw_boundserror(itp, (x,)))
    k = searchsortedfirst(itp.knots, x)
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
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        m[k] = Δk
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::FiniteDifferenceMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        else
            m[k] = (Δ[k-1] + Δk) / 2
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
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        else
            m[k] = (oneunit(TTension) - method.tension) * (y[k+1] - y[k-1]) / (x[k+1] - x[k-1])
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{TEl}, method::FritschCarlsonMonotonicInterpolation) where {TCoeffs, TEl}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        else
            # If any consecutive secant lines change sign (i.e. curve changes direction), initialize the tangent to zero.
            # This is needed to make the interpolation monotonic. Otherwise set tangent to the average of the secants.
            if Δk <= zero(Δk)
                m[k] = zero(TCoeffs)
            else
                m[k] = Δ[k-1] * (Δ[k-1] + Δk) / 2.0
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
    for k in 1:n-1
        Δk = Δ[k]
        if Δk == zero(TCoeffs)
            m[k] = zero(TCoeffs)
            m[k+1] = zero(TCoeffs)
            continue
        end
        α = m[k] / Δk
        β = m[k+1] / Δk
        τ = 3.0 / sqrt(α^2 + β^2)
        if τ < 1.0 # if we're outside the circle with radius 3 then move onto the circle
            m[k] = τ * α * Δk
            m[k+1] = τ * β * Δk
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
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        elseif Δ[k-1] * Δk <= zero(TCoeffs)
            m[k] = zero(TCoeffs)
        else
            α = (1.0 + (x[k+1] - x[k]) / (x[k+1] - x[k-1])) / 3.0
            m[k] = Δ[k-1] * Δk / (α*Δk + (1.0 - α)*Δ[k-1])
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
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        else
            p = ((x[k+1] - x[k]) * Δ[k-1] + (x[k] - x[k-1]) * Δk) / (x[k+1] - x[k-1])
            m[k] = (sign(Δ[k-1]) + sign(Δk)) *
                min(abs(Δ[k-1]), abs(Δk), 0.5*abs(p))
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{<:MonotonicInterpolationType}) = 1

lbounds(itp::MonotonicInterpolation) = (first(itp.knots),)
ubounds(itp::MonotonicInterpolation) = (last(itp.knots),)
