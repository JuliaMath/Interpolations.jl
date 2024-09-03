# specialized dispatch for 1D
function create_cubic_convolutional_coefs(knots::Knots1D, vs::AbstractVector{T}) where T
    N = length(vs) + 2
    c = zeros(T, N)
    c[2:N-1] = vs

    # 1D boundary conditions
    c[1] = 3*c[2] - 3*c[3] + c[4]
    c[N] = 3*c[N-1] - 3*c[N-2] + c[N-3]

    return c
end

# specialized dispatch for 2D
function create_cubic_convolutional_coefs(knots::Knots2D, vs::AbstractMatrix{T}) where T
    L, M = size(vs) .+ 2
    c = zeros(T, L, M)
    c[2:L-1, 2:M-1] = vs

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

# specialized dispatch for 3D
function create_cubic_convolutional_coefs(knots::Knots3D, vs::AbstractArray{T,3}) where {T}
    L, M, N = size(vs) .+ 2
    c = zeros(T, L, M, N)
    c[2:L-1, 2:M-1, 2:N-1] = vs

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

# generalized dispatch for N > 3
function create_cubic_convolutional_coefs(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N}) where {T,N}
    new_dims = size(vs) .+ 2    
    c = zeros(T, new_dims)
    inner_indices = map(d -> 2:(d-1), new_dims)
    c[inner_indices...] = vs
    
    # helper function to apply boundary condition along specific dimensions
    function apply_boundary_condition!(c, fixed_dims)
        for idx in CartesianIndices(size(c))
            is_boundary = any(idx[dim] in (1, size(c, dim)) for dim in fixed_dims)
            if !is_boundary
                continue
            end
            
            for dim in fixed_dims
                if idx[dim] == 1
                    offset1 = CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))
                    offset2 = CartesianIndex(ntuple(i -> i == dim ? 2 : 0, N))
                    offset3 = CartesianIndex(ntuple(i -> i == dim ? 3 : 0, N))
                    c[idx] = 3*c[idx + offset1] - 3*c[idx + offset2] + c[idx + offset3]
                elseif idx[dim] == size(c, dim)
                    offset1 = CartesianIndex(ntuple(i -> i == dim ? 1 : 0, N))
                    offset2 = CartesianIndex(ntuple(i -> i == dim ? 2 : 0, N))
                    offset3 = CartesianIndex(ntuple(i -> i == dim ? 3 : 0, N))
                    c[idx] = 3*c[idx - offset1] - 3*c[idx - offset2] + c[idx - offset3]
                end
            end
        end
    end
    
    # apply boundary conditions, starting from faces and working outwards towards corners
    for num_fixed_dims in 1:N
        for fixed_dims in combinations(1:N, num_fixed_dims)
            apply_boundary_condition!(c, fixed_dims)
        end
    end
    
    return c
end


# combinations function
function combinations(iterable, r)
    pool = collect(iterable)
    n = length(pool)
    if r > n
        return []
    end
    indices = collect(1:r)
    result = [pool[indices]]
    while true
        finished = true
        for i in reverse(1:r)
            if indices[i] != i + n - r
                finished = false
                indices[i] += 1
                for j in (i+1):r
                    indices[j] = indices[j-1] + 1
                end
                break
            end
        end
        if finished
            break
        end
        push!(result, pool[indices])
    end
    return result
end
