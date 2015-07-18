immutable NoInterp <: Degree{0} end

function define_indices_d(::Type{BSpline{NoInterp}}, d, pad)
    symix, symx = symbol("ix_",d), symbol("x_",d)
    :($symix = $symx)
end

function coefficients(::Type{BSpline{NoInterp}}, N, d)
    :()
end

function gradient_coefficients(::Type{BSpline{NoInterp}}, d)
    :()
end

function index_gen{IT<:DimSpec{BSpline}}(::Type{BSpline{NoInterp}}, ::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        return :($(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{BSpline{NoInterp}}, N) = 0
count_interp_dims{IT<:BSpline}(::Type{IT}, N) = N
function count_interp_dims{IT<:Tuple{Vararg{BSpline}}}(it::Type{IT}, N)
    n = 0
    for p in it.parameters
        n += count_interp_dims(p, 1)
    end
    n
end
