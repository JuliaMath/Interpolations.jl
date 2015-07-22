"""
`extrap_prep(T, Val{1}())`

1-dimensional, same lo/hi schemes
"""
extrap_prep{T}(::Type{T}, n::Val{1}) = extrap_prep(T, n, Val{1}())

"""
`extrap_prep(Tuple{T}, Val{1}())`

1-dimensional, same lo/hi schemes but specified as if N-dimensional

Equivalent with `extrap_prep(T, Val{1}())`
"""
extrap_prep{T}(::Type{Tuple{T}}, n::Val{1}) = extrap_prep(T, n)

"""
`extrap_prep(Tuple{T,T}, Val{1}())`

1-dimensional, same lo/hi schemes but specified as tuple

Equivalent with `extrap_prep(T, Val{1}())`
"""
extrap_prep{T}(::Type{Tuple{T,T}}, n::Val{1}) = extrap_prep(T, n)

"""
`extrap_prep(Tuple{Tuple{T,T}}, Val{1}())`

1-dimensional, same lo/hi schemes but specified as tuple and as if N-dimensional

Equivalent with `extrap_prep(T, Val{1}())`
"""
extrap_prep{T}(::Type{Tuple{Tuple{T,T}}}, n::Val{1}) = extrap_prep(T, n)

"""
`extrap_prep(Tuple{S,T}, Val{1}())`

1-dimensional, different lo/hi schemes
"""
function extrap_prep{S,T}(::Type{Tuple{S,T}}, n::Val{1})
    quote
        $(extrap_prep(S, n, Val{1}(), Val{:lo}()))
        $(extrap_prep(T, n, Val{1}(), Val{:hi}()))
    end
end

"""
`extrap_prep(Tuple{Tuple{S,T}}, Val{1}())`

1-dimensional, different lo/hi schemes but specified as if N-dimensional

Equivalent with `extrap_prep(Tuple{S,T}, Val{1}())`
"""
extrap_prep{S,T}(::Type{Tuple{Tuple{S,T}}}, n::Val{1}) = extrap_prep(Tuple{S,T}, n)

"""
`extrap_prep(Tuple{...}, Val{1}())`

1-dimensional, but incorrect tuple spec
""" # needed for ambiguity resolution
extrap_prep{T<:Tuple}(::Type{T}, ::Val{1}) = :(throw(ArgumentError("The 1-dimensional extrap configuration $T is not supported")))

"""
`extrap_prep(T, Val{N}())`

N-dimensional, all schemes same
"""
function extrap_prep{T,N}(::Type{T}, n::Val{N})
    exprs = Expr[]
    for d in 1:N
        push!(exprs, extrap_prep(T, n, Val{d}()))
    end
    return Expr(:block, exprs...)
end

"""
`extrap_prep(Tuple{...}, Val{N}())`

N-dimensional, different specs in different dimensions.

Note that the tuple must be of length N.
"""
function extrap_prep{N,T<:Tuple}(::Type{T}, n::Val{N})
    length(T.parameters) == N || return :(throw(ArgumentError("The $N-dimensional extrap configuration $T is not supported - must be a tuple of length $N (was length $(lenght(T.parameters)))")))
    exprs = Expr[]
    for d in 1:N
        Tdim = T.parameters[d]
        if Tdim <: Tuple
            length(Tdim.parameters) == 2 || return :(throw(ArgumentError("The extrap configuration $Tdim for dimension $d is not supported - must be a tuple of length 2 or a simple configuration type")))
            if Tdim.parameters[1] != Tdim.parameters[2]
                push!(exprs, extrap_prep(Tdim, n, Val{d}()))
            else
                push!(exprs, extrap_prep(Tdim.parameters[1], n, Val{d}()))
            end
        else
            push!(exprs, extrap_prep(Tdim, n, Val{d}()))
        end
    end
    return Expr(:block, exprs...)
end

"""
`extrap_prep(T, Val{N}(), Val{d}())`

N-dimensional fallback for expanding the same scheme lo/hi in a single dimension
"""
function extrap_prep{T,N,d}(::Type{T}, n::Val{N}, dim::Val{d})
    quote
        $(extrap_prep(T, n, dim, Val{:lo}()))
        $(extrap_prep(T, n, dim, Val{:hi}()))
    end
end

"""
`extrap_prep(Tuple{S,T}, Val{N}(), Val{d}())`

N-dimensional fallback for expanding the different schemes lo/hi in a single dimension
"""
function extrap_prep{S,T,N,d}(::Type{Tuple{S,T}}, n::Val{N}, dim::Val{d})
    quote
        $(extrap_prep(S, n, dim, Val{:lo}()))
        $(extrap_prep(T, n, dim, Val{:hi}()))
    end
end

"""
`extrap_prep(Tuple{T,T}, Val{N}(), Val{d}())`

N-dimensional fallback for expanding the same schemes lo/hi in a single dimension

Equivalent with `extrap_prep(T, Val{N}(), Val{d}())`
"""
function extrap_prep{T,N,d}(::Type{Tuple{T,T}}, n::Val{N}, dim::Val{d})
    quote
        $(extrap_prep(T, n, dim, Val{:lo}()))
        $(extrap_prep(T, n, dim, Val{:hi}()))
    end
end
