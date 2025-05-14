export CubicHermite

"""
    CubicHermite

This type is purposely left undocumented since the interface is expected to
radically change in order to make it conform to the `AbstractInterpolation`
interface. Consider this API to be highly unstable and non-public, and that
breaking changes to this code could be made in a point release.
"""
struct CubicHermite
    xs::Vector{Float64}
    ys::Vector{Float64}
    dydxs::Vector{Float64}
    function CubicHermite(xs, ys, dydxs)
        if length(xs) < 2
            throw(DomainError(
                "There must be at least two samples in the domain of the CubicHermite interpolator.",
            ))
        end
        x = xs[1]
        for i = 2:length(xs)
            if xs[i] ≤ x
                throw(DomainError("The list of abscissas must be strictly increasing."))
            end
            x = xs[i]
        end
        if (length(xs) != length(ys) || length(ys) != length(dydxs))
            throw(DomainError(
                "There must be exactly the same number of abscissas as ordinates and derivatives.",
            ))
        end
        new(xs, ys, dydxs)
    end
end

function (ch::CubicHermite)(x::Float64)::Float64
    if x ≤ ch.xs[1]
        if x == ch.xs[1]
            return ch.ys[1]
        end
        throw(DomainError(
            "Requested abscissas $x is less than the minimum value of the domain of the interpolator.",
        ))
    end

    if x ≥ last(ch.xs)
        if x == last(ch.xs)
            return last(ch.ys)
        end
        throw(DomainError(
            "Requested abscissa $x is greater than the maximum value of the domain of the interpolator.",
        ))
    end
    # searchsorted is convenient but it is a little too general for this usecase.
    # We have asserted that the xs's are strictly increasing.
    # Ostensibly there is another function that could do this more cleanly.
    urange = searchsorted(ch.xs, x)
    # urange.start exceeds urange.stop in this case which again indicates that
    # searchsorted is not quite the right tool for this job:
    i = urange.stop
    x1 = ch.xs[i]
    x2 = ch.xs[i+1]
    @assert x1 ≤ x ≤ x2
    h = x2 - x1
    y1 = ch.ys[i]
    y2 = ch.ys[i+1]

    dydx1 = ch.dydxs[i]
    dydx2 = ch.dydxs[i+1]
    # Corless, A Graduate Introduction to Numerical Methods, equation 14.10:
    c1 = (2x - 3x1 + x2) * (x - x2) * (x - x2) / (h * h * h)
    c2 = -(2x - 3x2 + x1) * (x - x1) * (x - x1) / (h * h * h)
    d1 = (x - x1) * (x - x2) * (x - x2) / (h * h)
    d2 = (x - x2) * (x - x1) * (x - x1) / (h * h)
    c1 * y1 + c2 * y2 + d1 * dydx1 + d2 * dydx2
end

function gradient(ch::CubicHermite, x::Float64)::Float64
    if x ≤ ch.xs[1]
        if x == ch.xs[1]
            return ch.dydxs[1]
        end
        throw(DomainError(
            "Requested abscissas $x is less than the minimum value of the domain of the interpolator.",
        ))
    end

    if x ≥ last(ch.xs)
        if x == last(ch.xs)
            return last(ch.dydxs)
        end
        throw(DomainError(
            "Requested abscissa $x is greater than the maximum value of the domain of the interpolator.",
        ))
    end
    urange = searchsorted(ch.xs, x)
    i = urange.stop
    x1 = ch.xs[i]
    x2 = ch.xs[i+1]
    @assert x1 ≤ x ≤ x2
    h = x2 - x1
    y1 = ch.ys[i]
    y2 = ch.ys[i+1]

    dydx1 = ch.dydxs[i]
    dydx2 = ch.dydxs[i+1]
    # Corless, A Graduate Introduction to Numerical Methods, equation 14.11:
    c1 = 6 * (x - x2) * (x - x1) / (h * h * h)
    d1 = (x - x2) * (3x - 2x1 - x2) / (h * h)
    d2 = (x - x1) * (3x - 2x2 - x1) / (h * h)
    c1 * (y1 - y2) + d1 * dydx1 + d2 * dydx2
end

function hessian(ch::CubicHermite, x::Float64)::Float64
    if x < ch.xs[1]
        throw(DomainError(
            "Requested abscissas $x is less than the minimum value of the domain of the interpolator.",
        ))
    end

    if x ≥ last(ch.xs)
        if x == last(ch.xs)
            h = ch.xs[end] - ch.xs[end-1]
            c1 = 6 / (h * h)
            d1 = 2 / h
            d2 = 4 / h
            return c1 * (ch.ys[end-1] - ch.ys[end]) +
                   d1 * ch.dydxs[end-1] +
                   d2 * ch.dydxs[end]
        end
        throw(DomainError(
            "Requested abscissa $x is greater than the maximum value of the domain of the interpolator.",
        ))
    end
    urange = searchsorted(ch.xs, x)
    i = urange.stop
    x1 = ch.xs[i]
    x2 = ch.xs[i+1]
    @assert x1 ≤ x ≤ x2
    h = x2 - x1
    y1 = ch.ys[i]
    y2 = ch.ys[i+1]

    dydx1 = ch.dydxs[i]
    dydx2 = ch.dydxs[i+1]
    # Corless, A Graduate Introduction to Numerical Methods, equation 14.12:
    c1 = 6 * (2x - x1 - x2) / (h * h * h)
    d1 = 2 * (3x - 2x2 - x1) / (h * h)
    d2 = 2 * (3x - 2x1 - x2) / (h * h)
    c1 * (y1 - y2) + d1 * dydx1 + d2 * dydx2
end
