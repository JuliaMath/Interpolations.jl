
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
    MonotonicInterpolationType,
    LinearMonotonicInterpolation,
    FiniteDifferenceMonotonicInterpolation,
    CardinalMonotonicInterpolation,
    FritschCarlsonMonotonicInterpolation,
    FritschButlandMonotonicInterpolation,
    SteffenMonotonicInterpolation

abstract type MonotonicInterpolationType <: InterpolationType end

struct LinearMonotonicInterpolation <: MonotonicInterpolationType
end

struct FiniteDifferenceMonotonicInterpolation <: MonotonicInterpolationType
end

struct CardinalMonotonicInterpolation{T<:Number} <: MonotonicInterpolationType
    tension :: T # must be in [0, 1]
end

struct FritschCarlsonMonotonicInterpolation <: MonotonicInterpolationType
end

struct FritschButlandMonotonicInterpolation <: MonotonicInterpolationType
end

struct SteffenMonotonicInterpolation <: MonotonicInterpolationType
end

struct MonotonicInterpolation{T, TCoeffs, Tel, Type<:MonotonicInterpolationType,
    K<:AbstractVector{<:Number}, AType <: AbstractArray{Tel,1}} <: AbstractInterpolation{T,1,DimSpec{Type}}

    it::Type
    knots::K
    A::AType
    m::Vector{TCoeffs}
    c::Vector{TCoeffs}
    d::Vector{TCoeffs}
end


size(A::MonotonicInterpolation) = size(A.knots)
axes(A::MonotonicInterpolation) = axes(A.knots)

function MonotonicInterpolation(::Type{TWeights}, it::IType, knots::K, A::AbstractArray{Tel,1},
    m::Vector{TCoeffs}, c::Vector{TCoeffs}, d::Vector{TCoeffs}) where {TWeights, TCoeffs, Tel, IType<:MonotonicInterpolationType, K<:AbstractVector{<:Number}}

    isconcretetype(IType) || error("The b-spline type must be a leaf type (was $IType)")
    isconcretetype(TCoeffs) || warn("For performance reasons, consider using an array of a concrete type (eltype(A) == $(eltype(A)))")

    check_monotonic(knots, A)

    cZero = zero(TWeights)
    if isempty(A)
        T = Base.promote_op(*, typeof(cZero), eltype(A))
    else
        T = typeof(cZero * first(A))
    end

    MonotonicInterpolation{T, TCoeffs, Tel, IType, K, typeof(A)}(it, knots, A, m, c, d)
end

function interpolate(::Type{TWeights}, ::Type{TCoeffs}, knots::K,
    A::AbstractArray{Tel,1}, it::IT) where {TWeights,TCoeffs,Tel,K<:AbstractVector{<:Number},IT<:MonotonicInterpolationType}

    check_monotonic(knots, A)

    # first we need to determine tangents (m)
    n = length(knots)
    m, Δ = calcTangents(TCoeffs, knots, A, it)
    c = Vector{TCoeffs}(undef, n-1)
    d = Vector{TCoeffs}(undef, n-1)
    for k ∈ 1:n-1
        if IT == LinearMonotonicInterpolation
            c[k] = d[k] = zero(TCoeffs)
        else
            xdiff = knots[k+1] - knots[k]
            c[k] = (3*Δ[k] - 2*m[k] - m[k+1]) / xdiff
            d[k] = (m[k] + m[k+1] - 2*Δ[k]) / (xdiff * xdiff)
        end
    end

    MonotonicInterpolation(TWeights, it, knots, A, m, c, d)
end

function interpolate(knots::AbstractVector{<:Number}, A::AbstractArray{Tel,1},
    it::IT) where {Tel,N,IT<:MonotonicInterpolationType}

    interpolate(tweight(A), tcoef(A), knots, A, it)
end

@inline function findKnot(itp::MonotonicInterpolation, x::Number)
    x >= itp.knots[1] || error("Given number $x is outside of interpolated range [$(itp.knots[1]), $(itp.knots[end])]")
    x <= itp.knots[end] || error("Given number $x is outside of interpolated range [$(itp.knots[1]), $(itp.knots[end])]")

    k = 1
    n = length(itp.knots)
    while k < n-1 && x > itp.knots[k+1]
        k += 1
    end
    return k
end

function (itp::MonotonicInterpolation)(x::Number)
    k = findKnot(itp, x)
    xdiff = x - itp.knots[k]
    return itp.A[k] + itp.m[k]*xdiff + itp.c[k]*xdiff*xdiff + itp.d[k]*xdiff*xdiff*xdiff
end

function derivative(itp::MonotonicInterpolation, x::Number)
    k = findKnot(itp, x)
    xdiff = x - itp.knots[k]
    return itp.m[k] + 2*itp.c[k]*xdiff + 3*itp.d[k]*xdiff*xdiff
end

@inline function check_monotonic(knots, A)
    axes(knots) == axes(A) || throw(DimensionMismatch("knot vector must have the same axes as the corresponding array"))
    issorted(knots) || error("knot-vector must be sorted in increasing order")
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{Tel}, method::LinearMonotonicInterpolation) where {TCoeffs, Tel}

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
    y::AbstractVector{Tel}, method::FiniteDifferenceMonotonicInterpolation) where {TCoeffs, Tel}

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
    y::AbstractVector{Tel}, method::CardinalMonotonicInterpolation{T}) where {T, TCoeffs, Tel}

    n = length(x)
    Δ = Vector{TCoeffs}(undef, n-1)
    m = Vector{TCoeffs}(undef, n)
    for k in 1:n-1
        Δk = (y[k+1] - y[k]) / (x[k+1] - x[k])
        Δ[k] = Δk
        if k == 1   # left endpoint
            m[k] = Δk
        else
            m[k] = (oneunit(T) - method.tension) * (y[k+1] - y[k-1]) / (x[k+1] - x[k-1])
        end
    end
    m[n] = Δ[n-1]
    return (m, Δ)
end

function calcTangents(::Type{TCoeffs}, x::AbstractVector{<:Number},
    y::AbstractVector{Tel}, method::FritschCarlsonMonotonicInterpolation) where {TCoeffs, Tel}

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
    y :: AbstractVector{Tel}, method :: FritschButlandMonotonicInterpolation) where {TCoeffs, Tel}

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
    y::AbstractVector{Tel}, method::SteffenMonotonicInterpolation) where {TCoeffs, Tel}

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

lbounds(itp::MonotonicInterpolation) = first(itp.knots)
ubounds(itp::MonotonicInterpolation) = last(itp.knots)
