import Base.LinAlg.LU, Base.getindex

### Tridiagonal inversion along a particular dimension, first offsetting the values by b

function A_ldiv_B_md!{T}(dest, F::LU{T,Tridiagonal{T}}, src, dim::Integer, b::AbstractVector)
    1 <= dim <= max(ndims(dest),ndims(src)) || throw(DimensionMismatch("The chosen dimension $dim is larger than $(ndims(src)) and $(ndims(dest))"))
    n = size(F, 1)
    n == size(src, dim) && n == size(dest, dim) || throw(DimensionMismatch("Sizes $n, $(size(src,dim)), and $(size(dest,dim)) do not match"))
    size(dest) == size(src) || throw(DimensionMismatch("Sizes $(size(dest)), $(size(src)) do not match"))
    check_matrix(F)
    length(b) == size(src,dim) || throw(DimensionMismatch("length(b) = $(length(b)), which does not match $(size(src,dim))"))
    R1 = CartesianRange(size(dest)[1:dim-1])  # these are not type-stable, so let's use a function barrier
    R2 = CartesianRange(size(dest)[dim+1:end])
    _A_ldiv_B_md!(dest, F, src, R1, R2, b)
end

# Filtering along the first dimension
function _A_ldiv_B_md!{T,CI<:CartesianIndex{0}}(dest, F::LU{T,Tridiagonal{T}}, src,  R1::CartesianRange{CI}, R2, b)
    dinv = 1./F.factors.d  # might not want to do this for small R2
    for I2 in R2
        invert_column!(dest, F, src, I2, b, dinv)
    end
    dest
end

function invert_column!{T}(dest, F::LU{T,Tridiagonal{T}}, src, I2, b, dinv)
    n = size(F, 1)
    if n == 0
        return nothing
    end
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    @inbounds begin
        # Forward substitution
        dest[1, I2] = src[1, I2] + b[1]
        for i = 2:n  # Can't use @simd here
            dest[i, I2] = src[i, I2] + b[i] - dl[i-1]*dest[i-1, I2]
        end
        # Backward substitution
        dest[n, I2] *= dinv[n]
        for i = n-1:-1:1   # Can't use @simd here
            dest[i, I2] = (dest[i, I2] - du[i]*dest[i+1, I2])*dinv[i]
        end
    end
    nothing
end

# Filtering along any other dimension
function _A_ldiv_B_md!{T}(dest, F::LU{T,Tridiagonal{T}}, src, R1, R2, b)
    n = size(F, 1)
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    dinv = 1./d  # might not want to do this for small R1 and R2
    # Forward substitution
    @inbounds for I2 in R2
        @simd for I1 in R1
            dest[I1, 1, I2] = src[I1, 1, I2] + b[1]
        end
        for i = 2:n
            @simd for I1 in R1
                dest[I1, i, I2] = src[I1, i, I2] + b[i] - dl[i-1]*dest[I1, i-1, I2]
            end
        end
        # Backward substitution
        @simd for I1 in R1
            dest[I1, n, I2] *= dinv[n]
        end
        for i = n-1:-1:1
            @simd for I1 in R1
                dest[I1, i, I2] = (dest[I1, i, I2] - du[i]*dest[I1, i+1, I2])*dinv[i]
            end
        end
    end
    dest
end

### Woodbury inversion along dimension 1
# Here we support only dimension 1 because the matrix multiplications in Woodbury
# inversion cannot be done in a cache-friendly way for any other dimension.
# It's therefore better to permutedims.
function filter_dim1!(dest, W::Woodbury, src, b)
    size(src,1) == size(W,1) == length(b) || throw(DimensionMismatch("Sizes $(size(src,1)), $(size(W,1)), and $(length(b)) do not match"))
    size(src) == size(dest) || throw(DimensionMismatch("Sizes $(size(dest)), $(size(src)) do not match"))
    check_matrix(W.A)
    R = CartesianRange(size(dest)[2:end])
    _filter_dim1!(dest, W, src, b, R)
end

function _filter_dim1!(dest, W::Woodbury, src, b, R)
    dinv = 1./W.A.factors.d
    n = length(dinv)
    tmp = W.tmpN2
    dl, du = W.A.factors.dl, W.A.factors.du
    for I in R   # iterate over "columns"
        invert_column!(dest, W.A, src, I, b, dinv)
        # TODO? Manually fuse the tridiagonal inversions with their "adjacent" matrix multiplications
        A_mul_b!(W.tmpk1, W.V, dest, I)
        A_mul_B!(W.tmpk2, W.Cp, W.tmpk1)
        A_mul_B!(tmp, W.U, W.tmpk2)
        # Fuses the final tridiagonal solve with the offset
        @inbounds begin
            # Forward substitution
            for i = 2:n  # Can't use @simd here
                tmp[i] = tmp[i] - dl[i-1]*tmp[i-1]
            end
            # Backward substitution
            t = tmp[n]*dinv[n]
            dest[n, I] -= t
            for i = n-1:-1:1   # Can't use @simd here
                t = (tmp[i] - du[i]*t)*dinv[i]
                dest[i, I] -= t
            end
        end
    end
    dest
end

function check_matrix{T}(F::LU{T,Tridiagonal{T}})
    for i = 1:size(F,1)
        F.ipiv[i] == i || error("For efficiency, pivoting is not supported")
    end
    du2 = F.factors.du2
    for i = 1:length(du2)
        du2[i] == 0 || error("For efficiency, du2 must be all zeros")
    end
end

check_matrix(M) = error("unsupported matrix type $(typeof(M))")

function A_mul_b!(dest, A, B, I)
    fill!(dest, 0)
    for j = 1:size(A,2)
        b = B[j,I]
        @simd for i = 1:size(A,1)
            dest[i] += A[i,j]*b
        end
    end 
end
