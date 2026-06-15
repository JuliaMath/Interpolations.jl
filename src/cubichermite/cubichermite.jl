#=
Assuming uniform knots with spacing 1, and taking the data and derivatives of
the function as given respectively by ({p_k}, {m_k}).

We can evaluate a point x \in (x_k, x_{k+1}) using the following formula

    p(x) = h00(x-xk) p_k + h10(x-xk) m_k + h01(x-xk) p_{k+1} + h11(x-xk) m_{k+1}

For more information, see:
    * https://en.wikipedia.org/wiki/Cubic_Hermite_spline
=#

"""
Build approximants to the tangent lines at the data points
"""
function approximate_derivatives(p)
    # Compute the deltas
    N = length(p)
    Δ = diff(p)

    # Initialize secants
    m = Array(Float64, N)
    m[1] = Δ[1]
    m[end] = Δ[N-1]
    for k in 2:N-1
        Δk = Δ[k]
        Δkm1 = Δ[k-1]
        # Average of secants unless opposite signs in which case use 0
        m[k] = Δk*Δkm1 > 0 ? (Δ[k-1] + Δ[k])/2 : 0.0
    end

    # Prevent overshoot and ensure monotonicity
    for k in 1:N-1
        # If the current point is the same as the next or previous point
        # Set derivative to 0
        if (abs(p[k+1] - p[k]) < 1e-14)
            m[k] = 0.0
            m[k+1] = 0.0

        # Otherwise, make sure that we satisfy the conditions for monotonicity
        # described in wiki page
        else
            # Compute α and β
            Δk = Δ[k]
            αk = m[k] / Δk
            βk = m[k+1] / Δk

            # Check whether sum of squares is less than 9
            if (αk^2 + βk^2) > 9
                # If it is, then fix it
                τ = 3 / sqrt(αk^2 + βk^2)
                m[k] = τ*αk*Δk
                m[k+1] = τ*βk*Δk
            end
            # More ad-hoc way to achieve the above
            # m[i] = αk < 3.0 ? αk : 3 * Δi
        end
    end

    return m
end

function pre_solve(p::Vector, m::Vector)
    N = length(p)
    @assert N == length(m)
    out = Array(Float64, 4, N-1)

    @inbounds @simd for n in 1:N-1
        pn = p[n]
        pn1 = p[n+1]
        mn = m[n]
        mn1 = m[n+1]
        out[1, n] = pn
        out[2, n] = mn
        out[3, n] = 3pn1 - 3pn - 2mn - mn1
        out[4, n] = -2pn1 + 2pn + mn + mn1
    end

    return out
end

function pre_solve(p::Vector)
    # Get number of points
    N = length(p)

    # Create approximate derivatives
    m = approximate_derivatives(p)
    c = pre_solve(p, m)

    return c
end



@inline function evaluate(coef_mat::Matrix, x::Real)
    k = floor(Int, x)
    t = x-k
    t2 = t*t
    t3 = t2*t

    return coef_mat[1, k] + coef_mat[2, k]*t + coef_mat[3, k]*t2 + coef_mat[4, k]*t3
end

