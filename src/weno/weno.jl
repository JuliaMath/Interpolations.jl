#=
    Interpolate using weighted essentially non-oscillatory (WENO) techniques.
=#

export WENO4Interpolation


"""
    WENOInterpolationType

Abstract class for all types of WENO interpolation.
"""
abstract type WENOInterpolationType <: InterpolationType end

"""
    WENO4Interpolation

Fourth-order WENO interpolation. Based on Janett et al. (2019),
"A novel fourth-order WENO interpolation technique. A possible new tool designed for radiative transfer",
https://ui.adsabs.harvard.edu/abs/2019A&A...624A.104J
"""
struct WENO4Interpolation{T} <: WENOInterpolationType where T <: Number
    nodes::Vector{T}
    α2::Vector{T}  # coefficients for the α2 and α3 weights
    α3::Vector{T}  # coefficients for the α2 and α3 weights
    q2::Array{T, 2}  # coefficients for the quadratic Lagrange polynomials q2
    q3::Array{T, 2}  # coefficients for the quadratic Lagrange polynomials q3
end

function Base.:(==)(o1::WENO4Interpolation, o2::WENO4Interpolation)
    o1.nodes == o2.nodes &&
    o1.α2 == o2.α2 &&
    o1.α3 == o2.α3 &&
    o1.q2 == o2.q2 &&
    o1.q3 == o2.q3
end

lbounds(itp::WENO4Interpolation) = (first(itp.nodes),)
ubounds(itp::WENO4Interpolation) = (last(itp.nodes),)


function interpolate(nodes::AbstractVector{<:Number}, A::AbstractVector{<:Number})
    n = length(nodes)
    α2 = Vector{eltype(nodes)}(undef, n)
    α3 = Vector{eltype(nodes)}(undef, n)
    q2 = Array{eltype(nodes)}(undef, 3, n)
    q3 = Array{eltype(nodes)}(undef, 3, n)

    for i ∈ 1:n
        xim, xi, xip, xipp = get_stencil_nodes(i, nodes)
        yim, yi, yip, yipp = get_stencil_values(i, A)
        q2[:, i], q3[:, i] = weno4_q(xim, xi, xip, xipp, yim, yi, yip, yipp)
        α2[i], α3[i] = weno4_α(xim, xi, xip, xipp, yim, yi, yip, yipp)
    end
    WENO4Interpolation(nodes, α2, α3, q2, q3)
end


function (itp::WENO4Interpolation)(x::Number)
    @boundscheck (checkbounds(Bool, itp, x) || Base.throw_boundserror(itp, (x,)))
    n = length(itp.nodes)
    i = binary_search(itp.nodes, x)
    xim, xi, xip, xipp = get_stencil_nodes(i, itp.nodes)
    Δim = x - xim
    Δi = x - xi
    Δip = x - xip
    Δipp = x - xipp
    q2 = itp.q2[1, i] * Δi * Δip + itp.q2[2, i] * Δim * Δip +  itp.q2[3, i] * Δim * Δi
    q3 = itp.q3[1, i] * Δip * Δipp + itp.q3[2, i] * Δi * Δipp +  itp.q3[3, i] * Δi * Δip
    if i == 1
        result = q3
    elseif i == n - 1
        result = q2
    else
        α2 = Δipp * itp.α2[i]
        α3 = Δip * itp.α3[i]
        ω2 = α2 / (α2 + α3)
        ω3 = α3 / (α2 + α3)
        result = ω2 * q2 + ω3 * q3
    end
    return result
end


function get_stencil_nodes(index::Int, nodes::AbstractVector{<:Number})
    n = length(nodes)
    if index == n
        index -= 1
    end
    xi = nodes[index]
    xip = nodes[index+1]
    if index == 1
        xim = zero(eltype(nodes))
        xipp = nodes[index+2]
    elseif index == n - 1
        xim = nodes[index-1]
        xipp = zero(eltype(nodes))
    else
        xim = nodes[index-1]
        xipp = nodes[index+2]
    end
    return (xim, xi, xip, xipp)
end


function get_stencil_values(index::Int, A::AbstractVector{<:Number})
    n = length(A)
    if index == n
        index -= 1
    end
    yi = A[index]
    yip = A[index+1]
    if index == 1
        yim = zero(eltype(A))
        yipp = A[index+2]
    elseif index == n - 1
        yim = A[index-1]
        yipp = zero(eltype(A))
    else
        yim= A[index-1]
        yipp = A[index+2]
    end
    return (yim, yi, yip, yipp)
end


"""
    binary_search(a::AbstractVector{<:Number}, x)

Return the index of the last value in `a` less than or equal to `x`.
`a` is assumed to be sorted.
"""
function binary_search(a::AbstractVector{<:Number}, x)
    n = length(a)
    ileft = 2
    i = n - 1
    while ileft <= i
        mid = (ileft + i) >> 1
        if @inbounds x < a[mid]
            i = mid - 1
        else
            ileft = mid + 1
        end
    end
    return i
end

"""
Compute the smoothness indicators β2 and β3 for 4th order WENO interpolation.
Arguments are the x and y values for the 4-point stencil.
"""
function weno4_α(xim, xi, xip, xipp, yim, yi, yip, yipp)
    ε = 1.0f-6
    him = xi - xim
    hi = xip - xi
    hip = xipp - xip
    H = him + hi + hip
    # Derivatives
    yyim = -((2 * him + hi) * H + him * (him + hi)) / (him * (him + hi) * H) * yim
    yyim += ((him + hi) * H) / (him * hi * (hi + hip)) * yi
    yyim -= (him * H) / ((him + hi) * hi * hip) * yip
    yyim += (him * (him + hi)) / ((hi + hip) * hip * H) * yipp
    yyi = -(hi * (hi + hip)) / (him * (him + hi) * H) * yim
    yyi += (hi * (hi + hip) - him * (2 * hi + hip)) / (him * hi * (hi + hip)) * yi
    yyi += (him * (hi + hip)) / ((him + hi) * hi * hip) * yip
    yyi -= (him * hi) / ((hi + hip) * hip * H) * yipp
    yyip = (hi * hip) / (him * (him + hi) * H) * yim
    yyip -= (hip * (him + hi)) / (him * hi * (hi + hip)) * yi
    yyip += ((him + 2 * hi) * hip - (him + hi) * hi) / ((him + hi) * hi * hip) * yip
    yyip += ((him + hi) * hi) / ((hi + hip) * hip * H) * yipp
    yyipp = -((hi + hip) * hip) / (him * (him + hi) * H) * yim
    yyipp += (hip * H) / (him * hi * (hi + hip)) * yi
    yyipp -= ((hi + hip) * H) / ((him + hi) * hi * hip) * yip
    yyipp += ((2 * hip + hi) * H + hip * (hi + hip)) / ((hi + hip) * hip * H) * yipp
    # Smoothness indicators
    β2 = (hi + hip)^2 * (abs(yyip - yyi) / hi - abs(yyi - yyim) / him)^2
    β3 = (him + hi)^2 * (abs(yyipp - yyip) / hip - abs(yyip - yyi) / hi)^2
    α2 = - one(typeof(xim)) / ((xipp - xim) * (ε + β2))
    α3 = one(typeof(xim)) / ((xipp - xim) * (ε + β3))
    return (α2, α3)
end

"""
Compute the q2 and q3 Lagrange interpolation factors for 4th order WENO interpolation.
Arguments are the x and y values for the 4-point stencil.
"""
function weno4_q(xim, xi, xip, xipp, yim, yi, yip, yipp)
    him = xi - xim
    hi = xip - xi
    hip = xipp - xip
    q2 = Vector{eltype(xim)}(undef, 3)
    q3 = Vector{eltype(xim)}(undef, 3)
    q2[1] = yim / (him * (him + hi))
    q2[2] = -yi / (him * hi)
    q2[3] = yip / ((him + hi) * hi)
    q3[1] = yi / (hi * (hi + hip))
    q3[2] = -yip / (hi * hip)
    q3[3] = yipp / ((hi + hip) * hip)
    return (q2, q3)
end
