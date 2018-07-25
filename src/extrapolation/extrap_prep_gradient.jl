# See ?extrap_prep for documentation for all these methods

extrap_prep(g::Val{:gradient}, ::Type{T}, n::Val{1}) where {T} = quote
    inds_etp = axes(etp)
    $(extrap_prep(g, T, n, Val{1}()))
end
extrap_prep(g::Val{:gradient}, ::Type{Tuple{T}}, n::Val{1}) where {T} = extrap_prep(g, T, n)
extrap_prep(g::Val{:gradient}, ::Type{Tuple{T,T}}, n::Val{1}) where {T} = extrap_prep(g, T, n)
extrap_prep(g::Val{:gradient}, ::Type{Tuple{Tuple{T,T}}}, n::Val{1}) where {T} = extrap_prep(g, T, n)
function extrap_prep(g::Val{:gradient}, ::Type{Tuple{S,T}}, ::Val{1}) where {S,T}
    quote
        inds_etp = axes(etp)
        $(extrap_prep(g, S, n, Val{1}(), Val{:lo}()))
        $(extrap_prep(g, T, n, Val{1}(), Val{:hi}()))
    end
end
extrap_prep(g::Val{:gradient}, ::Type{Tuple{Tuple{S,T}}}, n::Val{1}) where {S,T} = extrap_prep(g, Tuple{S,T}, n)
# needed for ambiguity resolution
extrap_prep(::Val{:gradient}, ::Type{T}, ::Val{1}) where {T<:Tuple} = :(throw(ArgumentError("The 1-dimensional extrap configuration $T is not supported")))


function extrap_prep(g::Val{:gradient}, ::Type{T}, n::Val{N}) where {T,N}
    Expr(:block, :(inds_etp = axes(etp)), [extrap_prep(g, T, n, Val{d}()) for d in 1:N]...)
end

function extrap_prep(g::Val{:gradient}, ::Type{T}, n::Val{N}) where {T<:Tuple,N}
    length(T.parameters) == N || return :(throw(ArgumentError("The $N-dimensional extrap configuration $T is not supported")))
    exprs = [:(inds_etp = axes(etp))]
    for d in 1:N
        Tdim = T.parameters[d]
        if Tdim <: Tuple
            length(Tdim.parameters) == 2 || return :(throw(ArgumentError("The extrap configuration $Tdim for dimension $d is not supported - must be a tuple of length 2 or a simple configuration type"))
                )
            if Tdim.parameters[1] != Tdim.parameters[2]
                push!(exprs, extrap_prep(g, Tdim, n, Val{d}()))
            else
                push!(exprs, extrap_prep(g, Tdim.parameters[1], n, Val{d}()))
            end
        else
            push!(exprs, extrap_prep(g, Tdim, n, Val{d}()))
        end
    end
    return Expr(:block, exprs...)
end

function extrap_prep(g::Val{:gradient}, ::Type{Tuple{S,T}}, n::Val{N}, dim::Val{d}) where {S,T,N,d}
    quote
        $(extrap_prep(g, S, n, dim, Val{:lo}()))
        $(extrap_prep(g, T, n, dim, Val{:hi}()))
    end
end

extrap_prep(g::Val{:gradient}, ::Type{T}, n::Val{N}, dim::Val{d}) where {T,N,d} = extrap_prep(g, Tuple{T,T}, n, dim)
extrap_prep(g::Val{:gradient}, args...) = extrap_prep(args...)
