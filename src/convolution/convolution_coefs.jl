function linear_predict(r)
    if length(r) == 1
        return [1.0] # constant
    elseif length(r) == 2
        return [2.0, -1.0] # linear
    elseif length(r) == 3
        return [3.0, -3.0, 1.0] # quadratic
    elseif length(r) == 4
        # third order boundary condition
        r0, r1, r2, r3 = r[1:4]
        coefs = zeros(3)
        coefs[1] = (r0^2*r1 - r0*r1*r2 - r0*r2*r3 - r1^3 + r1^2*r3 + r1*r2^2)/(isapprox((r0^3 - 2*r0*r1^2 - r0*r2^2 + 2*r1^2*r2), 0.0, atol=1e-6) ? 1e-6 : (r0^3 - 2*r0*r1^2 - r0*r2^2 + 2*r1^2*r2))
        coefs[2] = (r0*r2 - r1^2 - r1*r3 + r2^2)/(isapprox((r0^2 + r0*r2 - 2*r1^2), 0.0, atol=1e-6) ? 1e-6 : (r0^2 + r0*r2 - 2*r1^2))
        coefs[3] = (r0^2*r3 - 2*r0*r1*r2 + r1^3 - r1^2*r3 + r1*r2^2)/(isapprox((r0^3 - 2*r0*r1^2 - r0*r2^2 + 2*r1^2*r2), 0.0, atol=1e-6) ? 1e-6 : (r0^3 - 2*r0*r1^2 - r0*r2^2 + 2*r1^2*r2))
        return coefs
    elseif length(r) == 5
        # fourth order boundary condition
        r0, r1, r2, r3, r4 = r[1:5]
        coefs = zeros(4)
        coefs[1] = (r0^3*r1 - r0^2*r1*r2 - r0^2*r2*r3 - r0^2*r3*r4 - 2*r0*r1^3 + r0*r1^2*r3 + 2*r0*r1*r2*r4 + r0*r1*r3^2 + r0*r2^2*r3 + 3*r1^3*r2 - r1^3*r4 - 2*r1^2*r2*r3 + r1^2*r3*r4 - r1*r2^3 - r1*r2^2*r4 - r1*r2*r3^2 + r2^3*r3)/(r0^4 - 3*r0^2*r1^2 - 2*r0^2*r2^2 - r0^2*r3^2 + 4*r0*r1^2*r2 + 4*r0*r1*r2*r3 + r1^4 - 2*r1^3*r3 - 2*r1^2*r2^2 + r1^2*r3^2 - 2*r1*r2^2*r3 + r2^4)
        coefs[2] = (r0^3*r2 - r0^2*r1^2 - r0^2*r1*r3 - r0^2*r2*r4 + r0*r1^2*r4 + 3*r0*r1*r2*r3 + r0*r1*r3*r4 - r0*r2^3 - r0*r2*r3^2 + r1^4 - r1^3*r3 - r1^2*r2^2 - r1^2*r2*r4 - r1^2*r3^2 + 2*r1*r2^2*r3 - r1*r2*r3*r4 + r1*r3^3 + r2^3*r4 - r2^2*r3^2)/(r0^4 - 3*r0^2*r1^2 - 2*r0^2*r2^2 - r0^2*r3^2 + 4*r0*r1^2*r2 + 4*r0*r1*r2*r3 + r1^4 - 2*r1^3*r3 - 2*r1^2*r2^2 + r1^2*r3^2 - 2*r1*r2^2*r3 + r2^4)
        coefs[3] = (r0^3*r3 - 2*r0^2*r1*r2 - r0^2*r1*r4 + r0*r1^3 + 2*r0*r1*r2^2 + r0*r1*r2*r4 - r0*r2^2*r3 + r0*r2*r3*r4 - r0*r3^3 - r1^3*r2 + r1^3*r4 - 2*r1^2*r2*r3 - r1^2*r3*r4 + r1*r2^3 - r1*r2^2*r4 + 3*r1*r2*r3^2 - r2^3*r3)/(r0^4 - 3*r0^2*r1^2 - 2*r0^2*r2^2 - r0^2*r3^2 + 4*r0*r1^2*r2 + 4*r0*r1*r2*r3 + r1^4 - 2*r1^3*r3 - 2*r1^2*r2^2 + r1^2*r3^2 - 2*r1*r2^2*r3 + r2^4)
        coefs[4] = (r0^3*r4 - 2*r0^2*r1*r3 - r0^2*r2^2 + 3*r0*r1^2*r2 - 2*r0*r1^2*r4 + 2*r0*r1*r2*r3 - r0*r2^2*r4 + r0*r2*r3^2 - r1^4 + 2*r1^3*r3 - 2*r1^2*r2^2 + 2*r1^2*r2*r4 - r1^2*r3^2 - 2*r1*r2^2*r3 + r2^4)/(r0^4 - 3*r0^2*r1^2 - 2*r0^2*r2^2 - r0^2*r3^2 + 4*r0*r1^2*r2 + 4*r0*r1*r2*r3 + r1^4 - 2*r1^3*r3 - 2*r1^2*r2^2 + r1^2*r3^2 - 2*r1*r2^2*r3 + r2^4)
        return coefs
    else # if length(r) >= 6
        # fifth order boundary condition
        r0, r1, r2, r3, r4, r5 = r[1:6]
        coefs = zeros(5)
        coefs[1] = (r0^4*r1 - r0^3*r1*r2 - r0^3*r2*r3 - r0^3*r3*r4 - r0^3*r4*r5 - 3*r0^2*r1^3 + r0^2*r1^2*r3 - r0^2*r1*r2^2 + 2*r0^2*r1*r2*r4 + 2*r0^2*r1*r3*r5 + r0^2*r1*r4^2 + r0^2*r2^2*r3 + r0^2*r2^2*r5 + 2*r0^2*r2*r3*r4 + 6*r0*r1^3*r2 - r0*r1^3*r4 + 3*r0*r1^2*r2*r3 - 3*r0*r1^2*r2*r5 - r0*r1^2*r3*r4 + 2*r0*r1^2*r4*r5 - 4*r0*r1*r2^2*r4 - 4*r0*r1*r2*r3^2 - 2*r0*r1*r2*r3*r5 - r0*r1*r2*r4^2 - r0*r1*r3^2*r4 + r0*r2^2*r3*r4 + r0*r2^2*r4*r5 + r0*r2*r3^3 - r0*r2*r3^2*r5 - r0*r2*r3*r4^2 + r0*r3^3*r4 + r1^5 - 3*r1^4*r3 + r1^4*r5 - 5*r1^3*r2^2 + 2*r1^3*r2*r4 + 2*r1^3*r3^2 - 2*r1^3*r3*r5 - r1^3*r4^2 + r1^2*r2^2*r3 + 2*r1^2*r2^2*r5 + 2*r1^2*r2*r3*r4 - 2*r1^2*r2*r4*r5 + r1^2*r3^3 + r1^2*r3^2*r5 + r1^2*r3*r4^2 + 2*r1*r2^4 - r1*r2^2*r3^2 + 2*r1*r2^2*r3*r5 + r1*r2^2*r4^2 - 2*r1*r2*r3^2*r4 - r1*r3^4 - r2^4*r3 - r2^4*r5 + r2^2*r3^3)/(r0^5 - 4*r0^3*r1^2 - 3*r0^3*r2^2 - 2*r0^3*r3^2 - r0^3*r4^2 + 6*r0^2*r1^2*r2 + 8*r0^2*r1*r2*r3 + 4*r0^2*r1*r3*r4 + 2*r0^2*r2^2*r4 + 3*r0*r1^4 - 4*r0*r1^3*r3 - 2*r0*r1^2*r2^2 - 6*r0*r1^2*r2*r4 + 2*r0*r1^2*r4^2 - 8*r0*r1*r2^2*r3 - 4*r0*r1*r2*r3*r4 + 2*r0*r2^4 + 2*r0*r2^2*r3^2 + r0*r2^2*r4^2 - 2*r0*r2*r3^2*r4 + r0*r3^4 - 4*r1^4*r2 + 2*r1^4*r4 + 4*r1^3*r2*r3 - 4*r1^3*r3*r4 + 2*r1^2*r2^3 + 4*r1^2*r2^2*r4 + 4*r1^2*r2*r3^2 - 2*r1^2*r2*r4^2 + 2*r1^2*r3^2*r4 - 4*r1*r2^3*r3 + 4*r1*r2^2*r3*r4 - 4*r1*r2*r3^3 - 2*r2^4*r4 + 2*r2^3*r3^2)
        coefs[2] = (r0^4*r2 - r0^3*r1^2 - r0^3*r1*r3 - r0^3*r2*r4 - r0^3*r3*r5 - r0^2*r1^2*r2 + r0^2*r1^2*r4 + 3*r0^2*r1*r2*r3 + 2*r0^2*r1*r2*r5 + 3*r0^2*r1*r3*r4 + r0^2*r1*r4*r5 - 2*r0^2*r2^3 - r0^2*r2*r4^2 + 2*r0*r1^4 - r0*r1^3*r5 + 2*r0*r1^2*r2^2 - 4*r0*r1^2*r2*r4 - 3*r0*r1^2*r3^2 - r0*r1^2*r4^2 - 2*r0*r1*r2^2*r5 - r0*r1*r2*r3*r4 - r0*r1*r2*r4*r5 + r0*r1*r3^3 + r0*r1*r3*r4^2 + 4*r0*r2^3*r4 - r0*r2^2*r3^2 + r0*r2^2*r3*r5 - r0*r2*r3^2*r4 - r0*r2*r3*r4*r5 + r0*r2*r4^3 + r0*r3^3*r5 - r0*r3^2*r4^2 - 3*r1^4*r2 + r1^4*r4 + r1^3*r2*r3 + r1^3*r2*r5 - r1^3*r4*r5 + 2*r1^2*r2^3 + 2*r1^2*r2*r3^2 + 2*r1^2*r2*r3*r5 + 3*r1^2*r2*r4^2 - 2*r1^2*r3^2*r4 + r1^2*r3*r4*r5 - r1^2*r4^3 - 3*r1*r2^3*r3 - r1*r2^3*r5 - 3*r1*r2^2*r3*r4 + r1*r2^2*r4*r5 + r1*r2*r3^3 - 3*r1*r2*r3^2*r5 + r1*r2*r3*r4^2 + r1*r3^3*r4 + r2^3*r3^2 + r2^3*r3*r5 - 2*r2^3*r4^2 + 2*r2^2*r3^2*r4 - r2*r3^4)/(r0^5 - 4*r0^3*r1^2 - 3*r0^3*r2^2 - 2*r0^3*r3^2 - r0^3*r4^2 + 6*r0^2*r1^2*r2 + 8*r0^2*r1*r2*r3 + 4*r0^2*r1*r3*r4 + 2*r0^2*r2^2*r4 + 3*r0*r1^4 - 4*r0*r1^3*r3 - 2*r0*r1^2*r2^2 - 6*r0*r1^2*r2*r4 + 2*r0*r1^2*r4^2 - 8*r0*r1*r2^2*r3 - 4*r0*r1*r2*r3*r4 + 2*r0*r2^4 + 2*r0*r2^2*r3^2 + r0*r2^2*r4^2 - 2*r0*r2*r3^2*r4 + r0*r3^4 - 4*r1^4*r2 + 2*r1^4*r4 + 4*r1^3*r2*r3 - 4*r1^3*r3*r4 + 2*r1^2*r2^3 + 4*r1^2*r2^2*r4 + 4*r1^2*r2*r3^2 - 2*r1^2*r2*r4^2 + 2*r1^2*r3^2*r4 - 4*r1*r2^3*r3 + 4*r1*r2^2*r3*r4 - 4*r1*r2*r3^3 - 2*r2^4*r4 + 2*r2^3*r3^2)
        coefs[3] = (r0^2*r3 - 2*r0*r1*r2 - r0*r1*r4 + r0*r2*r3 - r0*r2*r5 + r0*r3*r4 + r1^3 + r1^2*r5 - 2*r1*r3^2 + r1*r3*r5 - r1*r4^2 + r2^2*r3 - r2^2*r5 + 2*r2*r3*r4 - r3^3)/(r0^3 + r0^2*r2 + r0^2*r4 - 3*r0*r1^2 - 2*r0*r1*r3 - 2*r0*r2^2 + r0*r2*r4 - r0*r3^2 + 4*r1^2*r2 - 2*r1^2*r4 + 4*r1*r2*r3 - 2*r2^3)
        coefs[4] = (r0^4*r4 - 2*r0^3*r1*r3 - r0^3*r1*r5 - r0^3*r2^2 + 3*r0^2*r1^2*r2 - r0^2*r1^2*r4 + 4*r0^2*r1*r2*r3 + r0^2*r1*r2*r5 - 2*r0^2*r2^2*r4 + r0^2*r2*r3^2 + r0^2*r2*r3*r5 - r0^2*r3^2*r4 + r0^2*r3*r4*r5 - r0^2*r4^3 - r0*r1^4 + r0*r1^3*r3 + 2*r0*r1^3*r5 - 4*r0*r1^2*r2^2 - 2*r0*r1^2*r3^2 - r0*r1^2*r3*r5 - 2*r0*r1*r2^2*r3 - 2*r0*r1*r2*r4*r5 + 2*r0*r1*r3^3 - r0*r1*r3^2*r5 + 3*r0*r1*r3*r4^2 + 2*r0*r2^4 - r0*r2^2*r3^2 - r0*r2^2*r3*r5 + 3*r0*r2^2*r4^2 - 2*r0*r2*r3^2*r4 + r1^4*r2 - r1^4*r4 + r1^3*r2*r3 - 3*r1^3*r2*r5 + 2*r1^3*r3*r4 + r1^3*r4*r5 + 4*r1^2*r2^2*r4 - 2*r1^2*r2*r3^2 + 2*r1^2*r2*r3*r5 - 3*r1^2*r2*r4^2 - 2*r1^2*r3^2*r4 - r1^2*r3*r4*r5 + r1^2*r4^3 - r1*r2^3*r3 + r1*r2^3*r5 + r1*r2^2*r3*r4 + r1*r2^2*r4*r5 + r1*r2*r3^3 + r1*r2*r3^2*r5 - 3*r1*r2*r3*r4^2 + r1*r3^3*r4 - 2*r2^4*r4 + r2^3*r3^2 - r2^3*r3*r5 + 2*r2^2*r3^2*r4 - r2*r3^4)/(r0^5 - 4*r0^3*r1^2 - 3*r0^3*r2^2 - 2*r0^3*r3^2 - r0^3*r4^2 + 6*r0^2*r1^2*r2 + 8*r0^2*r1*r2*r3 + 4*r0^2*r1*r3*r4 + 2*r0^2*r2^2*r4 + 3*r0*r1^4 - 4*r0*r1^3*r3 - 2*r0*r1^2*r2^2 - 6*r0*r1^2*r2*r4 + 2*r0*r1^2*r4^2 - 8*r0*r1*r2^2*r3 - 4*r0*r1*r2*r3*r4 + 2*r0*r2^4 + 2*r0*r2^2*r3^2 + r0*r2^2*r4^2 - 2*r0*r2*r3^2*r4 + r0*r3^4 - 4*r1^4*r2 + 2*r1^4*r4 + 4*r1^3*r2*r3 - 4*r1^3*r3*r4 + 2*r1^2*r2^3 + 4*r1^2*r2^2*r4 + 4*r1^2*r2*r3^2 - 2*r1^2*r2*r4^2 + 2*r1^2*r3^2*r4 - 4*r1*r2^3*r3 + 4*r1*r2^2*r3*r4 - 4*r1*r2*r3^3 - 2*r2^4*r4 + 2*r2^3*r3^2)
        coefs[5] = (r0^4*r5 - 2*r0^3*r1*r4 - 2*r0^3*r2*r3 + 3*r0^2*r1^2*r3 - 3*r0^2*r1^2*r5 + 3*r0^2*r1*r2^2 + 2*r0^2*r1*r2*r4 + r0^2*r1*r3^2 - 2*r0^2*r2^2*r5 + 2*r0^2*r2*r3*r4 - r0^2*r3^2*r5 + r0^2*r3*r4^2 - 4*r0*r1^3*r2 + 4*r0*r1^3*r4 - 2*r0*r1^2*r2*r3 + 4*r0*r1^2*r2*r5 - 2*r0*r1^2*r3*r4 - 2*r0*r1*r2^3 - 4*r0*r1*r2*r3^2 + 4*r0*r1*r2*r3*r5 - 2*r0*r1*r2*r4^2 - 2*r0*r1*r3^2*r4 + 2*r0*r2^3*r3 - 2*r0*r2^2*r3*r4 + 2*r0*r2*r3^3 + r1^5 - 3*r1^4*r3 + r1^4*r5 + 3*r1^3*r2^2 - 6*r1^3*r2*r4 + 2*r1^3*r3^2 - 2*r1^3*r3*r5 + r1^3*r4^2 + 5*r1^2*r2^2*r3 - 2*r1^2*r2^2*r5 + 4*r1^2*r2*r3*r4 + r1^2*r3^3 + r1^2*r3^2*r5 - r1^2*r3*r4^2 - 2*r1*r2^4 + 2*r1*r2^3*r4 - 5*r1*r2^2*r3^2 - 2*r1*r2^2*r3*r5 + r1*r2^2*r4^2 + 2*r1*r2*r3^2*r4 - r1*r3^4 + r2^4*r3 + r2^4*r5 - 2*r2^3*r3*r4 + r2^2*r3^3)/(r0^5 - 4*r0^3*r1^2 - 3*r0^3*r2^2 - 2*r0^3*r3^2 - r0^3*r4^2 + 6*r0^2*r1^2*r2 + 8*r0^2*r1*r2*r3 + 4*r0^2*r1*r3*r4 + 2*r0^2*r2^2*r4 + 3*r0*r1^4 - 4*r0*r1^3*r3 - 2*r0*r1^2*r2^2 - 6*r0*r1^2*r2*r4 + 2*r0*r1^2*r4^2 - 8*r0*r1*r2^2*r3 - 4*r0*r1*r2*r3*r4 + 2*r0*r2^4 + 2*r0*r2^2*r3^2 + r0*r2^2*r4^2 - 2*r0*r2*r3^2*r4 + r0*r3^4 - 4*r1^4*r2 + 2*r1^4*r4 + 4*r1^3*r2*r3 - 4*r1^3*r3*r4 + 2*r1^2*r2^3 + 4*r1^2*r2^2*r4 + 4*r1^2*r2*r3^2 - 2*r1^2*r2*r4^2 + 2*r1^2*r3^2*r4 - 4*r1*r2^3*r3 + 4*r1*r2^2*r3*r4 - 4*r1*r2*r3^3 - 2*r2^4*r4 + 2*r2^3*r3^2)
        return coefs
    end
end

function autocor_coefs(signal)
    acf = autocorrelation(signal)
    coefs = linear_predict(acf)
    if occursin("NaN", string(coefs))
        return [3.0, -3.0, 1.0] # fallback to quadratic
    end
    return coefs
end

# generalized dispatch all number of dimensions
function create_convolutional_coefs(vs::AbstractArray{T,N}, h::NTuple{N,T}, eqs::Int) where {T,N}
    new_dims = size(vs) .+ 2*(eqs-1)
    c = zeros(T, new_dims)
    inner_indices = map(d -> (1+(eqs-1)):(d-(eqs-1)), new_dims)
    c[inner_indices...] = vs

    # helper function to apply boundary condition along specific dimensions
    function apply_boundary_condition!(c, fixed_dims)
        for idx in CartesianIndices(size(c))
            is_boundary = any(idx[dim] in (1+(eqs - 1), size(c, dim) - (eqs - 1)) for dim in fixed_dims)
            if !is_boundary
                continue
            end

            for dim in fixed_dims

                if idx[dim] == 1 + (eqs - 1)

                    slice_offset = [CartesianIndex(ntuple(d -> d == dim ? i : 0, N)) for i = 0:size(vs, dim)-1]
                    slice = length(size(vs)) == 1 ? vs : [c[idx+slice_offset[i]] for i in 1:size(vs, dim)]
                    coef, y_offset, y_centered = boundary_coefs(slice, h[dim])
                    c_offset = [CartesianIndex(ntuple(d -> d == dim ? i : 0, N)) for i = 1:(eqs-1)]
                    y_extended = zeros((eqs-1)+length(y_centered))
                    y_extended[1+(eqs-1):end] .= y_centered
                    for j in 1:(eqs-1)
                        y_extended[1+(eqs-1)-j] = sum( coef[k] * y_extended[1 + (eqs-1) - j + k] for k = eachindex(coef))
                        c[idx-c_offset[j]] = y_offset + y_extended[1+(eqs-1)-j]
                    end

                elseif idx[dim] == size(c, dim) - (eqs - 1)

                    slice_offset = [CartesianIndex(ntuple(d -> d == dim ? i : 0, N)) for i = 0:size(vs, dim)-1]
                    slice = length(size(vs)) == 1 ? vs : [c[idx-slice_offset[i]] for i in 1:size(vs, dim)]
                    coef, y_offset, y_centered = boundary_coefs(slice, h[dim])
                    c_offset = [CartesianIndex(ntuple(d -> d == dim ? i : 0, N)) for i = 1:eqs-1]
                    y_extended = zeros(length(y_centered)+eqs-1)
                    y_extended[1:end-(eqs-1)] .= y_centered
                    for j in 1:eqs-1
                        y_extended[length(y_centered)+j] = sum( coef[k] * y_extended[length(y_centered) + j - k] for k = eachindex(coef))
                        c[idx+c_offset[j]] = y_offset + y_extended[length(y_centered)+j]
                    end    
                end
            end
        end
    end

    # Apply boundary conditions in order: 1 dimension, then 2, then 3, etc.
    for fixed_dims in 1:N
        apply_boundary_condition!(c, fixed_dims)
    end

    return c
end

function boundary_coefs(y, h::T) where T<:Number

    # check for linear signal
    n = length(y)
    y_mean = sum(y)/n
    x_uniform = 1:n
    x_mean = (n + 1) / 2
    x_centered = x_uniform .- x_mean
    y_centered = y .- y_mean    
    x_norm = sqrt(sum(x_centered.^2))
    y_norm = norm(y_centered)
    linear_score = abs(dot(x_centered, y_centered) / (x_norm * y_norm))
    if linear_score > 0.95
        return [2.0, -1.0, 0.0], y_mean, y_centered
    end

    # check for quadratic signal
    x_squared_centered = x_centered.^2 .- sum(x_centered.^2)/n
    quadratic_score = abs(dot(x_squared_centered, y_centered) / (norm(x_squared_centered) * y_norm))
    if quadratic_score > 0.95 && n > 5
        return [3.0, -3.0, 1.0], y_mean, y_centered
    end

    # assume periodic signal
    return autocor_coefs(y_centered), y_mean, y_centered

end

function autocorrelation(signal::Vector{T}) where T<:Number
    n = length(signal)
    signal_centered = signal .- sum(signal)/length(signal)
    variance = max(1e-6, sum(abs2, signal_centered) / n)
    acf = [sum(signal_centered[1:n-k] .* signal_centered[k+1:end]) for k in 0:n-1] ./ (n * variance)
    return acf
end