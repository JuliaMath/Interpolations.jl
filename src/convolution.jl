export CubicConvolutionalInterpolation, ConvolutionMethod

struct CubicConvolutionalKernel end

struct ConvolutionMethod <: InterpolationType end

function (::CubicConvolutionalKernel)(s)
    s_abs = abs(s)
    if s_abs < 1.0
        return (3/2)*s_abs^3 - (5/2)*s_abs^2 + 1
    elseif s_abs < 2.0
        return -(1/2)*s_abs^3 + (5/2)*s_abs^2 - 4*s_abs + 2
    else
        return 0.0
    end
end

struct CubicConvolutionalInterpolation{T,N,TCoefs<:AbstractArray,IT<:NTuple{N,ConvolutionMethod},Axs<:Tuple} <: AbstractInterpolation{T,N,IT}
    coefs::TCoefs
    knots::Axs
    it::IT
    h::NTuple{N,Float64}
    kernel::CubicConvolutionalKernel
end

function CubicConvolutionalInterpolation(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N}) where {T,N}
    if N > 3
        error("CubicConvolutionalInterpolation is currently implemented only for 1D, 2D and 3D. Got $N dimensions.")
    end

    coefs = create_cubic_convolutional_coefs(knots, vs)
    h = map(k -> k[2] - k[1], knots)
    it = ntuple(_ -> ConvolutionMethod(), N)
    if N == 1
        knots_new = (knots[1][1]-h[1]:h[1]:knots[1][end]+h[1],)
    elseif N == 2
        knots_new = (knots[1][1]-h[1]:h[1]:knots[1][end]+h[1], knots[2][1]-h[2]:h[2]:knots[2][end]+h[2])
    else
        knots_new = (knots[1][1]-h[1]:h[1]:knots[1][end]+h[1], knots[2][1]-h[2]:h[2]:knots[2][end]+h[2], knots[3][1]-h[3]:h[3]:knots[3][end]+h[3])
    end

    CubicConvolutionalInterpolation{T,N,typeof(coefs),typeof(it),typeof(knots)}(coefs, (knots_new), it, h, CubicConvolutionalKernel())
end

function create_cubic_convolutional_coefs(knots::Tuple{AbstractVector}, vs::AbstractVector{T}) where T
    N = length(vs) + 2
    c = zeros(T, N)
    c[2:N-1] = vs

    # 1D boundary conditions
    c[1] = 3*c[2] - 3*c[3] + c[4]
    c[N] = 3*c[N-1] - 3*c[N-2] + c[N-3]

    return c
end

function create_cubic_convolutional_coefs(knots::Tuple{AbstractVector,AbstractVector}, values::AbstractMatrix{T}) where T
    L, M = size(values) .+ 2
    c = zeros(T, L, M)
    c[2:L-1, 2:M-1] = values

    # 2D boundary conditions
    for i = 2:L-1
        c[i,1] = 3*c[i,2] - 3*c[i,3] + c[i,4]
        c[i,M] = 3*c[i,M-1] - 3*c[i,M-2] + c[i,M-3]
    end
    for j = 2:M-1
        c[1,j] = 3*c[2,j] - 3*c[3,j] + c[4,j]
        c[L,j] = 3*c[L-1,j] - 3*c[L-2,j] + c[L-3,j]
    end

    c[1,1] = 3*c[2,1] - 3*c[3,1] + c[4,1]
    c[L,1] = 3*c[L-1,1] - 3*c[L-2,1] + c[L-3,1]
    c[1,M] = 3*c[1,M-1] - 3*c[1,M-2] + c[1,M-3]
    c[L,M] = 3*c[L-1,M] - 3*c[L-2,M] + c[L-3,M]

    return c
end

function create_cubic_convolutional_coefs(knots::Tuple{AbstractVector,AbstractVector,AbstractVector}, values::AbstractArray{T,3}) where {T}
    L, M, N = size(values) .+ 2
    c = zeros(T, L, M, N)
    c[2:L-1, 2:M-1, 2:N-1] = values

    # 3D boundary conditions
    # 6 faces (1 fixed coordinate)
    for i = 2:L-1
        for j = 2:M-1
            c[i,j,1] = 3*c[i,j,2] - 3*c[i,j,3] + c[i,j,4]
            c[i,j,N] = 3*c[i,j,N-1] - 3*c[i,j,N-2] + c[i,j,N-3]
        end
    end
    for i = 2:L-1
        for k = 2:N-1
            c[i,1,k] = 3*c[i,2,k] - 3*c[i,3,k] + c[i,4,k]
            c[i,M,k] = 3*c[i,M-1,k] - 3*c[i,M-2,k] + c[i,M-3,k]
        end
    end
    for j = 2:M-1
        for k = 2:N-1
            c[1,j,k] = 3*c[2,j,k] - 3*c[3,j,k] + c[4,j,k]
            c[L,j,k] = 3*c[L-1,j,k] - 3*c[L-2,j,k] + c[L-3,j,k]
        end
    end
    # 12 edges (2 fixed coordinates)
    for i = 2:L-1
        c[i,1,1] = 3*c[i,2,1] - 3*c[i,3,1] + c[i,4,1]
        c[i,M,1] = 3*c[i,M-1,1] - 3*c[i,M-2,1] + c[i,M-3,1]
        c[i,1,N] = 3*c[i,1,N-1] - 3*c[i,1,N-2] + c[i,1,N-3]
        c[i,M,N] = 3*c[i,M,N-1] - 3*c[i,M,N-2] + c[i,M,N-3]
    end
    for j = 2:M-1
        c[1,j,1] = 3*c[2,j,1] - 3*c[3,j,1] + c[4,j,1]
        c[L,j,1] = 3*c[L-1,j,1] - 3*c[L-2,j,1] + c[L-3,j,1]
        c[1,j,N] = 3*c[1,j,N-1] - 3*c[1,j,N-2] + c[1,j,N-3]
        c[L,j,N] = 3*c[L,j,N-1] - 3*c[L,j,N-2] + c[L,j,N-3]
    end
    for k = 2:N-1
        c[1,1,k] = 3*c[2,1,k] - 3*c[3,1,k] + c[4,1,k]
        c[L,1,k] = 3*c[L-1,1,k] - 3*c[L-2,1,k] + c[L-3,1,k]
        c[1,M,k] = 3*c[1,M-1,k] - 3*c[1,M-2,k] + c[1,M-3,k]
        c[L,M,k] = 3*c[L,M-1,k] - 3*c[L,M-2,k] + c[L,M-3,k]
    end
    # 8 vertices (3 fixed coordinates)
    c[1,1,1] = 3*c[2,1,1] - 3*c[3,1,1] + c[4,1,1]
    c[L,1,1] = 3*c[L-1,1,1] - 3*c[L-2,1,1] + c[L-3,1,1]
    c[1,M,1] = 3*c[1,M-1,1] - 3*c[1,M-2,1] + c[1,M-3,1]
    c[L,M,1] = 3*c[L,M-1,1] - 3*c[L,M-2,1] + c[L,M-3,1]
    c[1,1,N] = 3*c[1,1,N-1] - 3*c[1,1,N-2] + c[1,1,N-3]
    c[L,1,N] = 3*c[L,1,N-1] - 3*c[L,1,N-2] + c[L,1,N-3]
    c[1,M,N] = 3*c[1,M,N-1] - 3*c[1,M,N-2] + c[1,M,N-3]
    c[L,M,N] = 3*c[L,M,N-1] - 3*c[L,M,N-2] + c[L,M,N-3]

    return c
end

function (itp::CubicConvolutionalInterpolation{T,1})(x::Number) where T
    i = searchsortedlast(itp.knots[1], x)

    i = limit_convolution_bounds(1, i, itp)

    s = (x-itp.knots[1][i])/itp.h[1]

    result = itp.coefs[i-1] * (-s^3 + 2*s^2 - s)/2 + itp.coefs[i] * (3*s^3 - 5*s^2 + 2)/2 +
            itp.coefs[i+1] * (-3*s^3 + 4*s^2 + s)/2 + itp.coefs[i+2] * (s^3 - s^2)/2  

    return result
end

function (itp::CubicConvolutionalInterpolation{T,2})(x::Number, y::Number) where T
    i = searchsortedlast(itp.knots[1], x)
    j = searchsortedlast(itp.knots[2], y)
    
    i = limit_convolution_bounds(1, i, itp)
    j = limit_convolution_bounds(2, j, itp)

    result = zero(T)
    for l = -1:2, m = -1:2
        result += itp.coefs[i+l, j+m] * 
                  itp.kernel((x - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((y - itp.knots[2][j+m]) / itp.h[2])
    end
    
    return result
end

function (itp::CubicConvolutionalInterpolation{T,3})(x::Number, y::Number, z::Number) where T
    i = searchsortedlast(itp.knots[1], x)
    j = searchsortedlast(itp.knots[2], y)
    k = searchsortedlast(itp.knots[3], z)

    i = limit_convolution_bounds(1, i, itp)
    j = limit_convolution_bounds(2, j, itp)
    k = limit_convolution_bounds(3, k, itp)

    result = zero(T)
    for l = -1:2, m = -1:2, n = -1:2
        result += itp.coefs[i+l, j+m, k+n] * 
                  itp.kernel((x - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((y - itp.knots[2][j+m]) / itp.h[2]) *
                  itp.kernel((z - itp.knots[3][k+n]) / itp.h[3])
    end
    
    return result
end

function limit_convolution_bounds(dim::Int, index::Int, itp::CubicConvolutionalInterpolation)
    if index < 2
        index = 2
    elseif index > length(itp.knots[dim]) - 2
        index = length(itp.knots[dim]) - 2
    end
    return index
end

Interpolations.getknots(itp::CubicConvolutionalInterpolation) = itp.knots
Base.axes(itp::CubicConvolutionalInterpolation) = axes(itp.coefs)
Base.size(itp::CubicConvolutionalInterpolation) = size(itp.coefs)
Interpolations.lbounds(itp::CubicConvolutionalInterpolation) = first.(itp.knots)
Interpolations.ubounds(itp::CubicConvolutionalInterpolation) = last.(itp.knots)
Interpolations.itpflag(::Type{<:CubicConvolutionalInterpolation{T,N,TCoefs,IT}}) where {T,N,TCoefs,IT} = IT()
Interpolations.coefficients(itp::CubicConvolutionalInterpolation) = itp.coefs
