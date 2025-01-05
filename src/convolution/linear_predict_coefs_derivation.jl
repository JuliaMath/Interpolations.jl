# small script used for symbolic derivation of autocorrelation coefficients
# this is done to avoid the linear solve when creating the boundary conditions

using SymPy

# function to set up the autocorrelation Toeplitz matrix
function toeplitz(c, r)
    n = length(r)
    T = zeros(eltype(c), n, n)
    for i in 1:n
        for j in 1:n
            if i >= j
                T[i,j] = c[i-j+1]
            else
                T[i,j] = r[j-i+1]
            end
        end
    end
    return T
end

function autocorrelation_matrix(p)
    # Create symbolic Toeplitz matrix
    r = SymPy.symbols("r:$(p)")
    R = toeplitz(r, r)
    
    return R
end

# Solve symbolically for the linear prediction coefficients for a given order
order = 4
R = autocorrelation_matrix(order)
r = R[2:end, 1]
linear_predict_coefs = R[1:end-1, 1:end-1] \ r