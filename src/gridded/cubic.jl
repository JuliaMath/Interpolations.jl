function getindex(itp::GriddedInterpolation{T,N,TCoefs,Interpolations.Gridded{Interpolations.Cubic{Interpolations.Natural}},K,P}, x::Number) where{T,N,TCoefs,K,P}
    a,b,c,d = coefficients(itp.knots[1], itp.coefs)
    interpolate(x, a, b, c, d, itp.knots[1])
end

function getindex(itp::GriddedInterpolation{T,N,TCoefs,Interpolations.Gridded{Interpolations.Cubic{Interpolations.Natural}},K,P}, x::AbstractVector) where{T,N,TCoefs,K,P}
    a,b,c,d = coefficients(itp.knots[1], itp.coefs)
    interpolate.(x, [a], [b], [c], [d], [itp.knots[1]])
end


"""
    interpolate(x, a, b, c, d, X, v=false)
Interpoalte at location x using coefficients a,b,c fitted to X & Y.
If i in in the range of x:
    Sᵢ(x) = a(x - xᵢ)³ + b(x - xᵢ)² + c(x - xᵢ) + dᵢ
Otherwise extrapolate with a = 0 (i.e. a 1st or 2nd order spline)
"""
function interpolate(x, a, b, c, d, X, v=false)
    idx = max(searchsortedfirst(X,x)-1,1)
    n   = length(X)
    idx = idx >= n ? idx - 1 : idx
    if idx < n                    # interpolation
        return a[idx]*(x^3) + b[idx]*(x^2) + c[idx]x + d[idx]
    end
    error()
end


## TODO: Check which one of the below is faster.
function expand_array!(x2t, x2)
    x2t[1]       = x2[1]
    x2t[end]     = x2[end]
    x2t[2:end-1] = cat(2, x2[2:end-1], x2[2:end-1])'[:]
    x2t
end

function expand_array2!(rhs, y)
    n = length(y)
    rhs[1] = y[1]        # f₀(1) = y[1] 
    for i in 2:(n-1)
        k = 2(i-1)       # fᵢ    = y[i]
        rhs[k]   = y[i]  # fᵢ₊₁  = y[i] 
        rhs[k+1] = y[i]
    end
    rhs[n] = y[n] # fₙ    = y[n]
end

function expand_array(x2::AbstractArray{T,1}) where T
    x2t = zeros(T,2(size(x2,1)-1))
    expand_array!(x2t, x2)
end

function coefficients(x, y;
        force_linear_extrapolation = true,
        boundary_condition = :natural)
    A,rhs = spline_coef_equations(x, y, force_linear_extrapolation=force_linear_extrapolation,
        boundary_condition=boundary_condition)
    coefficients = A\rhs
    a,b,c,d = [coefficients[n:4:end] for n in 1:4]
    return a,b,c,d
end

function spline_coef_equations(x, y;
        force_linear_extrapolation = true,
        boundary_condition = :natural
    )
    valid_boundary_conditions = [:natural, :periodic, :notaknot, :quadratic]
    if !(boundary_condition in valid_boundary_conditions)
        error("Boundary Condition must be one of $valid_boundary_conditions not $boundary_condition")
    end

    n   = length(x)
    A   = zeros(4(n-1),4(n-1))
    rhs = zeros(4(n-1))
    
    # First 2(n-1) polynomials are of the form:
    #     fᵢ = aᵢx³ + bᵢx² +cᵢx + dᵢ
    rhs[1:2(n-1)] = expand_array(y)

    x2     = expand_array(x)
    x2inds = expand_array(1:length(x))
    niind = 0
    for xi in 1:(length(x2))
        for ni in 1:4
            A[xi, (niind+ni)] = x2[xi]^(4-ni)
        end
        iseven(xi) && (niind += 4 )
    end

    # Next polynomials from first derivative
    #    3ax² + 2bx + c + 0
    niind = 1
    for xi in 2:(n-1)
        rind = xi + 2(n-1) - 1
        A[rind,niind:niind+3] .= A[rind,(niind+4):(niind+7)] .= [(3*(x[xi]^2)),2x[xi],1,0]
        A[rind,niind+4:niind+7] *= -1.
        niind += 4
    end
    
    # Next polynomials from 2nd derivative
    #    6ax + 2b + 0 + 0
    niind = 1
    for xi in 2:(n-1)
        rind = xi + 2(n-1) + n - 3
        A[rind,niind:niind+3] .= A[rind,(niind+4):(niind+7)] .= [6*(x[xi]),2,0,0]
        A[rind,niind+4:niind+7] *= -1.
        niind += 4
    end
    
    # Next boundary conditions:
    # Natural Spline Satisfies:
    # 6a₁x₁   + 2b₁ = 0
    A[end-1,1] = 6x[1]
    A[end-1,2] = 2.
    # 6aₙxₙ₊₁ + 2bₙ = 0
    A[end,4(n-2)+1]  = 6x[end]
    A[end,4(n-2)+2]  = 2.
    
    return A,rhs
end

#####

function define_indices_d(::Type{Gridded{Cubic{Reflect}}}, d, pad)
    symix, symixp, symx = Symbol("ix_",d), Symbol("ixp_",d), Symbol("x_",d)
    quote
        $symix = clamp($symix, 1, size(itp, $d)-1)
        $symixp = $symix + 1
    end
end

function coefficients(::Type{Gridded{Cubic{Reflect}}}, N, d)
    symix, symixp, symx = Symbol("ix_",d), Symbol("ixp_",d), Symbol("x_",d)
    sym, symp, symfx = Symbol("c_",d), Symbol("cp_",d), Symbol("fx_",d)
    symk, symkix = Symbol("k_",d), Symbol("kix_",d)
    quote
        $symkix = $symk[$symix]
        $symfx = ($symx - $symkix)/($symk[$symixp] - $symkix)
        $sym = 1 - $symfx
        $symp = $symfx
    end
end

function gradient_coefficients(::Type{Gridded{Cubic{Reflect}}}, d)
    sym, symp = Symbol("c_",d), Symbol("cp_",d)
    symk, symix = Symbol("k_",d), Symbol("ix_",d)
    symixp = Symbol("ixp_",d)
    quote
        $symp = 1/($symk[$symixp] - $symk[$symix])
        $sym = - $symp
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen(::Type{Gridded{Cubic{Reflect}}}, ::Type{IT}, N::Integer, offsets...) where IT<:DimSpec{Gridded}
    if length(offsets) < N
        d = length(offsets)+1
        sym = Symbol("c_", d)
        symp = Symbol("cp_", d)
        return :($sym * $(index_gen(IT, N, offsets..., 0)) + $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
