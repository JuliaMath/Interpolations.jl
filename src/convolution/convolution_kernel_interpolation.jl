# specialized dispatch for 1D
function (itp::ConvolutionInterpolation{T,1,TCoefs,IT,Axs,KA,Val{1},Val{DG},Val{EQ}})(x::Vararg{AbstractFloat,1}) where {T,TCoefs,IT,Axs,KA,DG,EQ}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)

    result = zero(T)
    for l = -(itp.eqs-1):itp.eqs
        result += itp.coefs[i+l] * 
                  itp.kernel( (x[1] - itp.knots[1][i+l]) / itp.h[1] )
    end
    
    return result
end

# specialized dispatch for 2D
function (itp::ConvolutionInterpolation{T,2,TCoefs,IT,Axs,KA,Val{2},Val{DG},Val{EQ}})(x::Vararg{AbstractFloat,2}) where {T,TCoefs,IT,Axs,KA,DG,EQ}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)

    result = zero(T)
    for l = -(itp.eqs-1):itp.eqs, 
        m = -(itp.eqs-1):itp.eqs
        result += itp.coefs[i+l, j+m] * 
                  itp.kernel((x[1] - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((x[2] - itp.knots[2][j+m]) / itp.h[2])
    end
    
    return result
end

# specialized dispatch for 3D
function (itp::ConvolutionInterpolation{T,3,TCoefs,IT,Axs,KA,Val{3},Val{DG},Val{EQ}})(x::Vararg{AbstractFloat,3}) where {T,TCoefs,IT,Axs,KA,DG,EQ}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)
    k = limit_convolution_bounds(3, searchsortedlast(itp.knots[3], x[3]), itp)

    result = zero(T)
    for l = -(itp.eqs-1):itp.eqs, 
        m = -(itp.eqs-1):itp.eqs, 
        n = -(itp.eqs-1):itp.eqs
        result += itp.coefs[i+l, j+m, k+n] * 
                  itp.kernel((x[1] - itp.knots[1][i+l]) / itp.h[1]) * 
                  itp.kernel((x[2] - itp.knots[2][j+m]) / itp.h[2]) *
                  itp.kernel((x[3] - itp.knots[3][k+n]) / itp.h[3])
    end
    
    return result
end

# generalized dispatch for N > 3
function (itp::ConvolutionInterpolation{T,N,TCoefs,IT,Axs,KA,HigherDimension{N},Val{DG},Val{EQ}})(x::Vararg{AbstractFloat,N}) where {T,N,TCoefs,IT,KA,Axs,DG,EQ}

    pos_ids = ntuple(d -> limit_convolution_bounds(d, ntuple(d -> searchsortedlast(itp.knots[d], x[d]), N)[d], itp), N)

    result = zero(T)
    for offsets in Iterators.product(ntuple(_ -> -(itp.eqs-1):itp.eqs, N)...)
        result += itp.coefs[(pos_ids .+ offsets)...] * 
                  prod(itp.kernel((x[d] - itp.knots[d][pos_ids[d] + offsets[d]]) / itp.h[d]) for d in 1:N)
    end
    
    return result
end