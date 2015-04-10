immutable Constant <: Degree{0} end

function define_indices{GT<:GridType}(::Type{BSpline{Constant,None}}, ::Type{GT}, N)
    :(@nexprs $N d->(ix_d = clamp(round(Int, real(x_d)), 1, size(itp, d))))
end

function coefficients{GT<:GridType}(::Type{BSpline{Constant,None}}, ::Type{GT}, N)
    :(@nexprs $N d->($(coefficients(BSpline{Constant,None}, GT, N, :d))))
end

function coefficients{GT<:GridType}(::Type{BSpline{Constant,None}}, ::Type{GT}, N, d)
    sym, symx = symbol(string("c_",d)), symbol(string("x_",d))
    :($sym = 1)
end

function gradient_coefficients{GT<:GridType}(::Type{BSpline{Constant,None}}, ::Type{GT}, N, d)
    sym, symx = symbol(string("c_",d)), symbol(string("x_",d))
    :($sym = 0)
end

function index_gen{IT<:BSpline{Constant}}(::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        d = length(offsets)+1
        sym = symbol("c_"*string(d))
        return :($sym * $(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
