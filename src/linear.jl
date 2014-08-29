
# Linear interpolation needs no BC
Interpolation{T,N,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Type{Linear}, ::Type{EB}) = Interpolation(A, Linear, BCnone, EB)

function body_gen(::Type{Linear}, N)
    quote
        # ix_d is the index in dimension d of the nearest node *before* the interpolation point
        # fx_d is a parameter in [0,1] such that x_d = ix_d + fx_d
        @nexprs $N d->(ix_d = ifloor(x_d); fx_d = x_d - convert(typeof(x_d), ix_d))
        # ixp_d is the index in dimension d of the nearest node *after* the interpolation point
        @nexprs $N d->(ixp_d = ix_d + 1)
        @inbounds ret = $(index_gen(Linear, N))
        ret
    end
end

function body_gen(::Type{Linear}, ::Type{ExtrapError}, N)
    quote
        # if x_d is outside the domain, throw a BoundsError
        @nexprs $N d->(1 <= x_d <= size(itp,d) || throw(BoundsError()))

        # else, use the general indexing scheme for Linear
        $(body_gen(Linear, N))
    end
end

function body_gen(::Type{Linear}, ::Type{ExtrapNaN}, N)
    quote
        # if x_d is outside the domain, return NaN
        @nexprs $N d->(1 <= x_d <= size(itp,d) || return(nan(eltype(itp))))

        # else, use the general indexing scheme for Linear
        $(body_gen(Linear, N))
    end
end

function body_gen(::Type{Linear}, ::Type{ExtrapConstant}, N)
    quote
        # if x_d is outside the domain, move it to the domain edge
        @nexprs $N d->(x_d = clamp(x_d, 1, size(itp,d)))

        # then, use the general indexing scheme for Linear
        $(body_gen(Linear, N))
    end
end

body_gen{BC<:BoundaryCondition,EB<:ExtrapolationBehavior}(::Type{Linear}, ::Type{BC}, ::Type{EB}, N) = body_gen(Linear, EB, N)

function index_gen(::Type{Linear}, N::Integer, offsets...)
    if length(offsets) < N
        sym = symbol("fx_"*string(length(offsets)+1))
        return :((one($sym)-$sym) * $(index_gen(Linear, N, offsets..., 0)) + $sym * $(index_gen(Linear, N, offsets..., 1)))
    else
        indices = [offsets[d] == 0 ? symbol("ix_"*string(d)) : symbol("ixp_"*string(d)) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

