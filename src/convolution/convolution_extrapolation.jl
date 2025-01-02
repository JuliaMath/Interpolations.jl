
struct ConvolutionExtrapolation{T,N,ITPT<:AbstractConvolutionInterpolation,ET<:BoundaryCondition} <: AbstractExtrapolation{T,N,ITPT,NTuple{N,ConvolutionMethod}}
    itp::ITPT
    et::ET
end

function extrapolate(itp::AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,KA,DT}, et::ET) where {T,N,TCoefs,IT,Axs,KA,DT,ET<:BoundaryCondition}
    ConvolutionExtrapolation{T,N,typeof(itp),ET}(itp, et)
end

getknots(etp::ConvolutionExtrapolation) = getknots(etp.itp)
Base.size(etp::ConvolutionExtrapolation) = size(etp.itp.coefs)
Base.axes(etp::ConvolutionExtrapolation) = axes(etp.itp.coefs)

function (etp::ConvolutionExtrapolation{T,N,ITPT,ET})(x::Vararg{AbstractFloat,N}) where {T,N,ITPT,ET}
    itp = etp.itp
    knots = getknots(itp)
    
    # Check if the point is within bounds
    within_bounds = true
    for d in 1:N
        if x[d] > knots[d][end-(itp.eqs-1)] || x[d] < knots[d][itp.eqs]
            within_bounds = false
            break
        end
    end
    
    if within_bounds
        return itp(x...)
    else
        return extrapolate_point(etp, x)
    end
end

function extrapolate_point(etp::ConvolutionExtrapolation{T,N,ITPT,ET}, x::NTuple{N,AbstractFloat}) where {T,N,ITPT,ET}
    itp = etp.itp
    knots = getknots(itp)
    
    function reflect(y, l, u)
        yr = mod(y - l, 2(u-l)) + l
        return ifelse(yr > u, 2u-yr, yr)
    end

    clamped_x = zeros(T, N)
    for d in 1:N
        if x[d] < knots[d][itp.eqs]
            clamped_x[d] = knots[d][itp.eqs]
        elseif x[d] > knots[d][end-(itp.eqs-1)]
            clamped_x[d] = knots[d][end-(itp.eqs-1)]
        else
            clamped_x[d] = x[d]
        end
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
                period = knots[i][end-(itp.eqs-1)] - knots[i][itp.eqs]
                x_periodic = knots[i][itp.eqs] + mod(x[i] - knots[i][itp.eqs], period)
                val = itp(ntuple(j -> j == i ? x_periodic : clamped_x[j], N)...)
            elseif etp.et isa Reflect
                x_reflect = reflect(x[i], knots[i][itp.eqs], knots[i][end-(itp.eqs-1)])
                val = itp(ntuple(j -> j == i ? x_reflect : clamped_x[j], N)...)
            elseif etp.et isa Throw
                continue
            else
                throw(BoundsError(itp, x))
            end
        end
    end
    
    return val
end

function calculate_gradient(etp::ConvolutionExtrapolation{T,N,ITPT,ET}, x::NTuple{N,AbstractFloat}) where {T,N,ITPT,ET}
    itp = etp.itp
    knots = getknots(itp)
    grad = zeros(T, N)

    for d in 1:N
        grad[d] = calculate_component_gradient(etp, x, d)
    end

    return grad
end

function calculate_component_gradient(etp::ConvolutionExtrapolation{T,N,ITPT,ET}, x::NTuple{N,AbstractFloat}, d::Int) where {T,N,ITPT,ET}
    itp = etp.itp
    knots = getknots(itp)
    h = etp.itp.h[d] * 1e-6  # Small step size

    # Helper function to evaluate interpolator along dimension d
    function eval_along_dim(offset)
        return itp(ntuple(i -> i == d ? x[i] + offset : x[i], N)...)
    end

    return (eval_along_dim(2h) - 8 * eval_along_dim(h) + 8 * eval_along_dim(-h) - eval_along_dim(-2h)) / (12h)

end

bounds(etp::ConvolutionExtrapolation) = bounds(etp.itp)