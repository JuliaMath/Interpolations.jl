struct CubicConvolutionalExtrapolation{T,N,ITPT<:CubicConvolutionalInterpolation,ET<:BoundaryCondition,O} <: AbstractExtrapolation{T,N,ITPT,NTuple{N,ConvolutionMethod}}
    itp::ITPT
    et::ET
end

function Interpolations.extrapolate(itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,KA,DT,O}, et::ET) where {T,N,TCoefs,IT,Axs,KA,DT,O,ET<:BoundaryCondition}
    CubicConvolutionalExtrapolation{T,N,typeof(itp),ET,O}(itp, et)
end

Interpolations.getknots(etp::CubicConvolutionalExtrapolation) = getknots(etp.itp)
Base.size(etp::CubicConvolutionalExtrapolation) = size(etp.itp.coefs)
Base.axes(etp::CubicConvolutionalExtrapolation) = axes(etp.itp.coefs)

function (etp::CubicConvolutionalExtrapolation{T,N,ITPT,ET,O})(x::Vararg{Number,N}) where {T,N,ITPT,ET,O}
    itp = etp.itp
    knots = getknots(itp)
    
    # Check if the point is within bounds
    within_bounds = true
    if O == Interpolations.OrderOfAccuracy{3}
        for d in 1:N
            if x[d] > knots[d][end-1] || x[d] < knots[d][2]
                within_bounds = false
                break
            end
        end
    elseif O == Interpolations.OrderOfAccuracy{4}
        for d in 1:N
            if x[d] > knots[d][end-2] || x[d] < knots[d][3]
                within_bounds = false
                break
            end
        end
    else
        throw(ArgumentError("Order of accuracy $O not supported"))
    end
    
    if within_bounds
        return itp(x...)
    else
        return extrapolate_point(etp, x)
    end
end

function extrapolate_point(etp::CubicConvolutionalExtrapolation{T,N,ITPT,ET,O}, x::NTuple{N,Number}) where {T,N,ITPT,ET,O}
    itp = etp.itp
    knots = getknots(itp)
    
    clamped_x = zeros(T, N)
    if O == Interpolations.OrderOfAccuracy{3}
        for d in 1:N
            if x[d] < knots[d][2]
                clamped_x[d] = knots[d][2]
            elseif x[d] > knots[d][end-1]
                clamped_x[d] = knots[d][end-1]
            else
                clamped_x[d] = x[d]
            end
        end
    elseif O == Interpolations.OrderOfAccuracy{4}
        for d in 1:N
            if x[d] < knots[d][3]
                clamped_x[d] = knots[d][3]
            elseif x[d] > knots[d][end-2]
                clamped_x[d] = knots[d][end-2]
            else
                clamped_x[d] = x[d]
            end
        end
    else
        throw(ArgumentError("Order of accuracy $O not supported"))
    end
    clamped_x = ntuple(i -> clamped_x[i], N)
    grad = calculate_gradient(etp, clamped_x)
    val = itp(clamped_x...)

    for i in 1:N
        if !isapprox(x[i], clamped_x[i])
            if etp.et isa Line
                val += grad[i] * (x[i] - clamped_x[i])
            elseif etp.et isa Flat
                # Do nothing, keep the boundary value
            elseif etp.et isa Periodic
                if O == Interpolations.OrderOfAccuracy{3}
                    period = knots[i][end-1] - knots[i][2]
                    x_periodic = knots[i][2] + mod(x[i] - knots[i][2], period)    
                elseif O == Interpolations.OrderOfAccuracy{4}
                    period = knots[i][end-2] - knots[i][3]
                    x_periodic = knots[i][3] + mod(x[i] - knots[i][3], period)    
                else
                    throw(ArgumentError("Order of accuracy $O not supported"))
                end
                val = itp(ntuple(j -> j == i ? x_periodic : clamped_x[j], N)...)
            elseif etp.et isa Reflect
                if O == Interpolations.OrderOfAccuracy{3}
                    x_reflect = reflect(x[i], knots[i][2], knots[i][end-1])
                elseif O == Interpolations.OrderOfAccuracy{4}
                    x_reflect = reflect(x[i], knots[i][3], last(knots[i][end-2]))
                else
                    throw(ArgumentError("Order of accuracy $O not supported"))
                end
                val = itp(ntuple(j -> j == i ? x_reflect : clamped_x[j], N)...)
            elseif etp.et isa Throw
                throw(BoundsError(itp, x))
            end
        end
    end
    
    return val
end

function calculate_gradient(etp::CubicConvolutionalExtrapolation{T,N,ITPT,ET,O}, x::NTuple{N,Number}) where {T,N,ITPT,ET,O}
    itp = etp.itp
    knots = getknots(itp)
    grad = zeros(T, N)

    for d in 1:N
        grad[d] = calculate_component_gradient(etp, x, d)
    end

    return grad
end

function calculate_component_gradient(etp::CubicConvolutionalExtrapolation{T,N,ITPT,ET,O}, x::NTuple{N,Number}, d::Int) where {T,N,ITPT,ET,O}
    itp = etp.itp
    knots = getknots(itp)
    h = etp.itp.h[d] * 1e-6  # Small step size

    # Helper function to evaluate interpolator along dimension d
    function eval_along_dim(offset)
        return itp(ntuple(i -> i == d ? x[i] + offset : x[i], N)...)
    end

    if O == Interpolations.OrderOfAccuracy{3}
        if x[d] ≈ knots[d][2]  # Left boundary
            return (-11 * eval_along_dim(0) + 18 * eval_along_dim(h) - 9 * eval_along_dim(2h) + 2 * eval_along_dim(3h)) / (6h)
        elseif x[d] ≈ knots[d][end-1]  # Right boundary
            return (11 * eval_along_dim(0) - 18 * eval_along_dim(-h) + 9 * eval_along_dim(-2h) - 2 * eval_along_dim(-3h)) / (6h)
        else  # Interior point
            return (eval_along_dim(2h) - 8 * eval_along_dim(h) + 8 * eval_along_dim(-h) - eval_along_dim(-2h)) / (12h)
        end
    elseif O == Interpolations.OrderOfAccuracy{4}
        if x[d] ≈ knots[d][3]  # Left boundary
            return (-11 * eval_along_dim(0) + 18 * eval_along_dim(h) - 9 * eval_along_dim(2h) + 2 * eval_along_dim(3h)) / (6h)
        elseif x[d] ≈ knots[d][end-2]  # Right boundary
            return (11 * eval_along_dim(0) - 18 * eval_along_dim(-h) + 9 * eval_along_dim(-2h) - 2 * eval_along_dim(-3h)) / (6h)
        else  # Interior point
            return (eval_along_dim(2h) - 8 * eval_along_dim(h) + 8 * eval_along_dim(-h) - eval_along_dim(-2h)) / (12h)
        end
    end
end

function reflect(y, l, u)
    yr = mod(y - l, 2(u-l)) + l
    return ifelse(yr > u, 2u-yr, yr)
end

Interpolations.bounds(etp::CubicConvolutionalExtrapolation) = bounds(etp.itp)