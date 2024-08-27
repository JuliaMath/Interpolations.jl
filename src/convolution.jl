export ConvolutionInterpolation, ConvolutionMethod

struct ConvolutionKernel end

struct ConvolutionMethod <: InterpolationType end

function (::ConvolutionKernel)(s)
    s_abs = abs(s)
    if s_abs < 1.0
        return (3/2)*s_abs^3 - (5/2)*s_abs^2 + 1
    elseif s_abs < 2.0
        return -(1/2)*s_abs^3 + (5/2)*s_abs^2 - 4*s_abs + 2
    else
        return 0.0
    end
end

struct ConvolutionInterpolation{T,N,TCoefs<:AbstractArray,IT<:NTuple{N,ConvolutionMethod},Axs<:Tuple} <: AbstractInterpolation{T,N,IT}
    coefs::TCoefs
    knots::Axs
    it::IT
    h::NTuple{N,Float64}
    kernel::ConvolutionKernel
end

function ConvolutionInterpolation(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N}) where {T,N}
    if N > 2
        error("ConvolutionInterpolation is currently implemented only for 1D and 2D. Got $N dimensions.")
    end

    coefs = create_coefs(knots, vs)
    h = map(k -> k[2] - k[1], knots)
    it = ntuple(_ -> ConvolutionMethod(), N)
    if N == 1
        knots_new = (knots[1][1]-h[1]:h[1]:knots[1][end]+h[1],)
    else
        knots_new = (knots[1][1]-h[1]:h[1]:knots[1][end]+h[1], knots[2][1]-h[2]:h[2]:knots[2][end]+h[2])
    end

    ConvolutionInterpolation{T,N,typeof(coefs),typeof(it),typeof(knots)}(coefs, (knots_new), it, h, ConvolutionKernel())
end

function create_coefs(knots::Tuple{AbstractVector}, vs::AbstractVector{T}) where T
    N = length(vs) + 2
    c = zeros(T, N)
    c[2:N-1] = vs

    # 1D boundary conditions
    c[1] = 3*c[2] - 3*c[3] + c[4]
    c[N] = 3*c[N-1] - 3*c[N-2] + c[N-3]

    return c
end

function create_coefs(knots::Tuple{AbstractVector,AbstractVector}, values::AbstractMatrix{T}) where T
    N, M = size(values) .+ 2
    c = zeros(T, N, M)
    c[2:N-1, 2:M-1] = values

    # 2D boundary conditions
    for k = 2:M-1
        c[1,k] = 3*c[2,k] - 3*c[3,k] + c[4,k]
        c[N,k] = 3*c[N-1,k] - 3*c[N-2,k] + c[N-3,k]
    end
    for j = 2:N-1
        c[j,1] = 3*c[j,2] - 3*c[j,3] + c[j,4]
        c[j,M] = 3*c[j,M-1] - 3*c[j,M-2] + c[j,M-3]
    end
    c[1,1] = 3*c[2,1] - 3*c[3,1] + c[4,1]
    c[N,1] = 3*c[N-1,1] - 3*c[N-2,1] + c[N-3,1]
    c[1,M] = 3*c[1,M-1] - 3*c[1,M-2] + c[1,M-3]
    c[N,M] = 3*c[N-1,M] - 3*c[N-2,M] + c[N-3,M]

    return c
end

function (itp::ConvolutionInterpolation{T,1})(x::Number) where T
    j = searchsortedlast(itp.knots[1], x)

    if j == 1
        j = j + 1
    elseif j == length(itp.knots[1])-1
        j = j - 1
    elseif j == length(itp.knots[1])
        j = j - 2
    end

    s = (x-itp.knots[1][j])/itp.h[1]

    result = itp.coefs[j-1] * (-s^3 + 2*s^2 - s)/2 + itp.coefs[j] * (3*s^3 - 5*s^2 + 2)/2 +
            itp.coefs[j+1] * (-3*s^3 + 4*s^2 + s)/2 + itp.coefs[j+2] * (s^3 - s^2)/2  

    return result
end

function (itp::ConvolutionInterpolation{T,2})(x::Number, y::Number) where T
    j = searchsortedlast(itp.knots[1], x)
    k = searchsortedlast(itp.knots[2], y)
    
    if j < 1 || j == 1
        j = 2
    elseif j == length(itp.knots[1])-1
        j = j - 1
    elseif j == length(itp.knots[1])
        j = j - 2
    end

    if k < 1 || k == 1
        k = 2
    elseif k == length(itp.knots[2])-1
        k = k - 1
    elseif k == length(itp.knots[2])
        k = k - 2
    end

    result = zero(T)
    for l = -1:2, m = -1:2
        result += itp.coefs[j+l, k+m] * 
                  itp.kernel((x - itp.knots[1][j+l]) / itp.h[1]) * 
                  itp.kernel((y - itp.knots[2][k+m]) / itp.h[2])
    end
    
    return result
end

Interpolations.getknots(itp::ConvolutionInterpolation) = itp.knots
Base.axes(itp::ConvolutionInterpolation) = axes(itp.coefs)
Base.size(itp::ConvolutionInterpolation) = size(itp.coefs)
Interpolations.lbounds(itp::ConvolutionInterpolation) = first.(itp.knots)
Interpolations.ubounds(itp::ConvolutionInterpolation) = last.(itp.knots)
Interpolations.itpflag(::Type{<:ConvolutionInterpolation{T,N,TCoefs,IT}}) where {T,N,TCoefs,IT} = IT()
Interpolations.coefficients(itp::ConvolutionInterpolation) = itp.coefs

function Base.getindex(itp::ConvolutionInterpolation{T,1}, x::Number) where T
    itp(x)
end

function Base.getindex(itp::ConvolutionInterpolation{T,2}, x::Number, y::Number) where T
    itp(x, y)
end