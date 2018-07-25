function interpolate(A::AbstractArray, ::NoInterp, gt::GT) where {GT<:DimSpec{GridType}}
    interpolate(Int, eltype(A), A, NoInterp(), gt)
end

iextract(::Type{NoInterp}, d) = NoInterp

function define_indices_d(::Type{NoInterp}, d, pad)
    symix, symx = Symbol("ix_",d), Symbol("x_",d)
    :($symix = convert(Int, $symx))
end

function coefficients(::Type{NoInterp}, N, d)
    :()
end

function index_gen(::Type{NoInterp}, ::Type{IT}, N::Integer, offsets...) where IT<:DimSpec
    if (length(offsets) < N)
        return :($(index_gen(IT, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

# FIXME: needed due to a Julia inference bug (in julia 0.7)
function index_gen(::Type{NoInterp}, ::Type{NoInterp}, N::Integer, offsets...)
    if (length(offsets) < N)
        return :($(index_gen(NoInterp, NoInterp, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

padding(::Type{NoInterp}) = Val{0}()

# How many non-NoInterp dimensions are there?
count_interp_dims(::Type{NoInterp}, N) = 0
count_interp_dims(::Type{IT}, N) where {IT<:InterpolationType} = N
function count_interp_dims(it::Type{IT}, N) where IT<:Tuple{Vararg{InterpolationType}}
    n = 0
    for p in it.parameters
        n += count_interp_dims(p, 1)
    end
    n
end

prefilter(::Type{TWeights}, ::Type{TC}, A, ::Type{IT},::Type{GT}) where {TWeights, TC, IT<:NoInterp, GT<:GridType} = A, Val{0}()
