function define_indices_d(::Type{NoInterp}, d, pad)
    symix, symx = Symbol("ix_",d), Symbol("x_",d)
    :($symix = convert(Int, $symx))
end

function coefficients(::Type{NoInterp}, N, d)
    :()
end

function index_gen{IT<:DimSpec}(::Type{NoInterp}, ::Type{IT}, N::Integer, offsets...)
    if (length(offsets) < N)
        return :($(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

padding(::Type{NoInterp}) = Val{0}()

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{NoInterp}, N) = 0
count_interp_dims{IT<:InterpolationType}(::Type{IT}, N) = N
function count_interp_dims{IT<:Tuple{Vararg{InterpolationType}}}(it::Type{IT}, N)
    n = 0
    for p in it.parameters
        n += count_interp_dims(p, 1)
    end
    n
end
