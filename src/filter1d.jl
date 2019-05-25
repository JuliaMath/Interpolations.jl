import Base: getindex
import AxisAlgorithms: A_ldiv_B_md!, _A_ldiv_B_md!

### Tridiagonal inversion along a particular dimension, first offsetting the values by b

A_ldiv_B_md!(dest, ::Nothing, src, dim::Integer, ::Nothing) = dest

function A_ldiv_B_md!(dest, F, src, dim::Integer, b::AbstractVector)
    1 <= dim <= max(ndims(dest),ndims(src)) || throw(DimensionMismatch("The chosen dimension $dim is larger than $(ndims(src)) and $(ndims(dest))"))
    n = size(F, 1)
    n == size(src, dim) && n == size(dest, dim) || throw(DimensionMismatch("Sizes $n, $(size(src,dim)), and $(size(dest,dim)) do not match"))
    size(dest) == size(src) || throw(DimensionMismatch("Sizes $(size(dest)), $(size(src)) do not match"))
    check_matrix(F)
    length(b) == size(src,dim) || throw(DimensionMismatch("length(b) = $(length(b)), which does not match $(size(src,dim))"))
    R1 = CartesianIndices(size(dest)[1:dim-1])  # these are not type-stable, so let's use a function barrier
    R2 = CartesianIndices(size(dest)[dim+1:end])
    _A_ldiv_B_md!(dest, F, src, R1, R2, b)
end

# Filtering along the first dimension
function _A_ldiv_B_md!(dest, F::LU{T,<:Tridiagonal{T}}, src,  R1::CartesianIndices{0}, R2, b) where {T}
    n = size(F, 1)
    if n == 0
        return nothing
    end
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    dinv = 1 ./ d  # might not want to do this for small R2
    @inbounds for I2 in R2
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
function _A_ldiv_B_md!(dest, F::LU{T,<:Tridiagonal{T}}, src, R1, R2, b) where T
    n = size(F, 1)
    dl = F.factors.dl
    d  = F.factors.d
    du = F.factors.du
    dinv = 1 ./ d  # might not want to do this for small R1 and R2
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

### Woodbury solver with offset

function _A_ldiv_B_md!(dest, W::Woodbury, src,  R1, R2, b)
    _A_ldiv_B_md!(dest, W.A, src, R1, R2, b)
    tmp1 = AxisAlgorithms._A_mul_B_md(W.V, dest, R1, R2)
    tmp2 = AxisAlgorithms._A_mul_B_md(W.Cp, tmp1, R1, R2)
    tmp3 = AxisAlgorithms._A_mul_B_md(W.U, tmp2, R1, R2)
    # TODO?: would be nice to fuse the next two steps
    AxisAlgorithms._A_ldiv_B_md!(tmp3, W.A, tmp3, R1, R2)
    AxisAlgorithms.sub!(dest, tmp3)
end

function check_matrix(F::LU{T,<:Tridiagonal{T}}) where T
    for i = 1:size(F,1)
        F.ipiv[i] == i || error("For efficiency, pivoting is not supported")
    end
    du2 = F.factors.du2
    for i = 1:length(du2)
        du2[i] == 0 || error("For efficiency, du2 must be all zeros")
    end
end

check_matrix(F::Woodbury) = check_matrix(F.A)

check_matrix(M) = error("unsupported matrix type $(typeof(M))")
