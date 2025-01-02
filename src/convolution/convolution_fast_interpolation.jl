# specialized dispatch for 1D
function (itp::FastConvolutionInterpolation{T,1,TCoefs,IT,Axs,KA,Val{1},Val{DG},EQ,PR,KP})(x::Vararg{AbstractFloat,1}) where {T,TCoefs,IT,Axs,KA,DG,EQ,PR,KP}

    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    x_diff_left = (x[1] - itp.knots[1][i])/itp.h[1] # Normalized distance to left sample
    x_diff_right = 1 - x_diff_left # Normalized distance to right sample
    index_nearest = searchsortednearest(itp.pre_range, x_diff_right) # lookup index in precomputed matrix
    result = sum(itp.coefs[i+j]*itp.kernel_pre[index_nearest, j+itp.eqs] for j in -(itp.eqs-1):itp.eqs)
    return result
end

# specialized dispatch for 2D
function (itp::FastConvolutionInterpolation{T,2,TCoefs,IT,Axs,KA,Val{2},Val{DG},EQ,PR,KP})(x::Vararg{AbstractFloat,2}) where {T,TCoefs,IT,Axs,KA,DG,EQ,PR,KP}

    # first dimension
    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    x_diff_left = (x[1] - itp.knots[1][i])/itp.h[1] # Normalized distance to left sample
    x_diff_right = 1 - x_diff_left # Normalized distance to right sample
    index_nearest_x = searchsortednearest(itp.pre_range, x_diff_right) # lookup index in precomputed matrix
    # second dimension
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)
    y_diff_left = (x[2] - itp.knots[2][j])/itp.h[2] # Normalized distance to left sample
    y_diff_right = 1 - y_diff_left # Normalized distance to right sample
    index_nearest_y = searchsortednearest(itp.pre_range, y_diff_right) # lookup index in precomputed matrix

    # result
    result = sum(sum(itp.coefs[i+l, j+m]*itp.kernel_pre[index_nearest_x, l+itp.eqs]*itp.kernel_pre[index_nearest_y, m+itp.eqs] for l in -(itp.eqs-1):itp.eqs) for m in -(itp.eqs-1):itp.eqs)

    return result
end

# specialized dispatch for 3D
function (itp::FastConvolutionInterpolation{T,3,TCoefs,IT,Axs,KA,Val{3},Val{DG},EQ,PR,KP})(x::Vararg{AbstractFloat,3}) where {T,TCoefs,IT,Axs,KA,DG,EQ,PR,KP}

    # first dimension
    i = limit_convolution_bounds(1, searchsortedlast(itp.knots[1], x[1]), itp)
    x_diff_left = (x[1] - itp.knots[1][i])/itp.h[1] # Normalized distance to left sample
    x_diff_right = 1 - x_diff_left # Normalized distance to right sample
    index_nearest_x = searchsortednearest(itp.pre_range, x_diff_right) # lookup index in precomputed matrix
    # second dimension
    j = limit_convolution_bounds(2, searchsortedlast(itp.knots[2], x[2]), itp)
    y_diff_left = (x[2] - itp.knots[2][j])/itp.h[2] # Normalized distance to left sample
    y_diff_right = 1 - y_diff_left # Normalized distance to right sample
    index_nearest_y = searchsortednearest(itp.pre_range, y_diff_right) # lookup index in precomputed matrix
    # third dimension
    k = limit_convolution_bounds(3, searchsortedlast(itp.knots[3], x[3]), itp)
    z_diff_left = (x[3] - itp.knots[3][k])/itp.h[3] # Normalized distance to left sample
    z_diff_right = 1 - z_diff_left # Normalized distance to right sample
    index_nearest_z = searchsortednearest(itp.pre_range, z_diff_right) # lookup index in precomputed matrix

    # result
    result = sum(sum(sum(itp.coefs[i+l, j+m, k+n]*itp.kernel_pre[index_nearest_x, l+itp.eqs]*
                                                itp.kernel_pre[index_nearest_y, m+itp.eqs]*
                                                itp.kernel_pre[index_nearest_z, n+itp.eqs]
                                                for l in -(itp.eqs-1):itp.eqs) for m in -(itp.eqs-1):itp.eqs) for n in -(itp.eqs-1):itp.eqs )

    return result
end

# generalized dispatch for N > 3
function (itp::FastConvolutionInterpolation{T,N,TCoefs,IT,Axs,KA,HigherDimension{N},Val{DG},EQ,PR,KP})(x::Vararg{AbstractFloat,N}) where {T,N,TCoefs,IT,KA,Axs,DG,EQ,PR,KP}

    pos_ids = ntuple(d -> limit_convolution_bounds(d, ntuple(d -> searchsortedlast(itp.knots[d], x[d]), N)[d], itp), N)
    diff_left = ntuple(d -> (x[d] - itp.knots[d][pos_ids[d]])/itp.h[d], N) # Normalized distances to left sample
    diff_right = ntuple(d -> 1 - diff_left[d], N) # Normalized distance to right sample
    index_nearest = ntuple(d -> searchsortednearest(itp.pre_range, diff_right[d]), N) # lookup index in precomputed matrix

    result = zero(T)
    for offsets in Iterators.product(ntuple(_ -> -(itp.eqs-1):itp.eqs, N)...)
        result += itp.coefs[(pos_ids .+ offsets)...] * prod(itp.kernel_pre[index_nearest[d], offsets[d]+itp.eqs] )
    end
    
    return result
end

###
function limit_convolution_bounds(dim::Int, index::Int, itp::AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},Val{N},Val{DG},EQ}) where {T,N,TCoefs,IT,Axs,DG,EQ}
    if index < itp.eqs
        index = itp.eqs
    elseif index > length(itp.knots[dim]) - itp.eqs
        index = length(itp.knots[dim]) - itp.eqs
    end
    return index
end

function limit_convolution_bounds(dim::Int, index::Int, itp::AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},HigherDimension{N},Val{DG},EQ}) where {T,N,TCoefs,IT,Axs,DG,EQ}
    if index < itp.eqs
        index = itp.eqs
    elseif index > length(itp.knots[dim]) - itp.eqs
        index = length(itp.knots[dim]) - itp.eqs
    end
    return index
end


function limit_convolution_bounds(dim::Int, index::Int, itp::AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,GaussianConvolutionKernel{B},Val{N},Val{DG},EQ}) where {T,N,TCoefs,IT,Axs,DG,EQ,B}
    if index < itp.eqs
        index = itp.eqs
    elseif index > length(itp.knots[dim]) - itp.eqs
        index = length(itp.knots[dim]) - itp.eqs
    end
    return index
end

function limit_convolution_bounds(dim::Int, index::Int, itp::AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,GaussianConvolutionKernel{B},HigherDimension{N},Val{DG},EQ}) where {T,N,TCoefs,IT,Axs,DG,EQ,B}
    if index < itp.eqs
        index = itp.eqs
    elseif index > length(itp.knots[dim]) - itp.eqs
        index = length(itp.knots[dim]) - itp.eqs
    end
    return index
end

function searchsortednearest(a,x)
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
        return idx
    else
        return idx-1
    end
end