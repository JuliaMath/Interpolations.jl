# specialized dispatch for 1D
function (itp::CubicConvolutionalInterpolation{T,1,TCoefs,IT,Axs,KA,Val{1},OrderOfAccuracy{O}})(x::Vararg{T,1}) where {T,TCoefs,IT,Axs,KA,O}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)

    result = zero(T)
    for l = (O == 3 ? (-1:2) : (-2:3))
        result += itp.coefs[i+l] * 
                  itp.kernel((x[1] - itp.knots[1][i+l]) / itp.h[1])
    end
    
    return result
end

# specialized dispatch for 2D
function (itp::CubicConvolutionalInterpolation{T,2,TCoefs,IT,Axs,KA,Val{2},OrderOfAccuracy{O}})(x::Vararg{T,2}) where {T,TCoefs,IT,Axs,KA,O}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)

    result = zero(T)
    for l = (O == 3 ? (-1:2) : (-2:3)), m = (O == 3 ? (-1:2) : (-2:3))
        result += itp.coefs[i+l, j+m] * 
                  itp.kernel((x[1] - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((x[2] - itp.knots[2][j+m]) / itp.h[2])
    end
    
    return result
end

# specialized dispatch for 3D
function (itp::CubicConvolutionalInterpolation{T,3,TCoefs,IT,Axs,KA,Val{3},OrderOfAccuracy{O}})(x::Vararg{T,3}) where {T,TCoefs,IT,Axs,KA,O}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)
    k = limit_convolution_bounds(3, searchsortedlast(itp.knots[3], x[3]), itp)

    result = zero(T)
    for l = (O == 3 ? (-1:2) : (-2:3)), m = (O == 3 ? (-1:2) : (-2:3)), n = (O == 3 ? (-1:2) : (-2:3))
        result += itp.coefs[i+l, j+m, k+n] * 
                  itp.kernel((x[1] - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((x[2] - itp.knots[2][j+m]) / itp.h[2]) *
                  itp.kernel((x[3] - itp.knots[3][k+n]) / itp.h[3])
    end
    
    return result
end

# generalized dispatch for N > 3
function (itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,KA,HigherDimension{N},OrderOfAccuracy{O}})(x::Vararg{T,N}) where {T,N,TCoefs,IT,KA,Axs,O}

    pos_ids = ntuple(d -> limit_convolution_bounds(d, ntuple(d -> searchsortedlast(itp.knots[d], x[d]), N)[d], itp), N)

    result = zero(T)
    for offsets in Iterators.product(ntuple(_ -> (O == 3 ? (-1:2) : (-2:3)), N)...)
        result += itp.coefs[(pos_ids .+ offsets)...] * 
                  prod(itp.kernel((x[d] - itp.knots[d][pos_ids[d] + offsets[d]]) / itp.h[d]) for d in 1:N)
    end
    
    return result
end

function limit_convolution_bounds(dim::Int, index::Int, itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{3},Val{N},OrderOfAccuracy{3}}) where {T,N,TCoefs,IT,Axs}
    if index < 2
        index = 2
    elseif index > length(itp.knots[dim]) - 2
        index = length(itp.knots[dim]) - 2
    end
    return index
end

function limit_convolution_bounds(dim::Int, index::Int, itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{3},HigherDimension{N},OrderOfAccuracy{3}}) where {T,N,TCoefs,IT,Axs}
    if index < 2
        index = 2
    elseif index > length(itp.knots[dim]) - 2
        index = length(itp.knots[dim]) - 2
    end
    return index
end

function limit_convolution_bounds(dim::Int, index::Int, itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{4},Val{N},OrderOfAccuracy{4}}) where {T,N,TCoefs,IT,Axs}
    if index < 3
        index = 3
    elseif index > length(itp.knots[dim]) - 3
        index = length(itp.knots[dim]) - 3
    end
    return index
end

function limit_convolution_bounds(dim::Int, index::Int, itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{4},HigherDimension{N},OrderOfAccuracy{4}}) where {T,N,TCoefs,IT,Axs}
    if index < 3
        index = 3
    elseif index > length(itp.knots[dim]) - 3
        index = length(itp.knots[dim]) - 3
    end
    return index
end
