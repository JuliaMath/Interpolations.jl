"""
The `extrap_prep` function is used by `getindex_impl` to generate the body of
the `getindex` function for extrapolation objects.

The methods of `extrap_prep` work in "layers", iteratively working out the exact
expression needed.

The first layer takes a specification of the extrapolation scheme(s) to be used
and a `Val` object that specifies the dimensionality of the extrapolation object:
`extrap_prep{T,N}(::Type{T}, Val{N})`. These methods only dispatch to the second
layer, and need not be extended for new schemes.

The second layer also takes a `Val` object that specifies a single dimension on
which to work: `extrap_prep{T,N,d}(::Type{T}, ::Val{N}, ::Val{d}). The methods
with this signature in src/extrapolation/extrap_prep.jl simply expand into a
block with sub-expressions for handling too-low and too-high values separately
(the third layer), but specific interpolation schemes can provide more specific
methods for this layer that handle both ends simultaneously. For example, the
`Flat` scheme has a layer-2 method that uses `clamp` to restrict the coordinate
when used in both directions, but uses `min` and `max` when handling each end
separately.

The third layer, to which the second dispatches if no scheme-specific method is
found, adds a final `Val` object with a symbol `:lo` or `:hi`:
`extrap_prep{T,N,d,l}(::Type{T}, ::Val{N}, ::Val{d}, ::Val{l})`. These methods
must be specified for each extrapolation scheme. However, the general framework
takes care of expanding all possible tuple combinations, so individual schemes
need only care about e.g. `T==Flat`.

In addition to these methods, there is a similar three-layer method hierarchy
for gradient evaluation, in which a `Val{:gradient}` is prepended to the other
arguments:
`extrap_prep{T,N,d,l}(::Val{:gradient}`, ::Type{T}, ::Val{N}, ::Val{d}, ::Val{l})`
If nothing else is specified for the individual schemes, these methods forward
to the same methods without the `:gradient` argument, i.e. the same behavior as
for value extrapolation. This works well with all schemes that are simple
coordinate transformations, but for anything else methods for the low- and high-
value cases need to be implemented for each scheme.
""" extrap_prep

extrap_prep{T}(::Type{T}, n::Val{1}) = extrap_prep(T, n, Val{1}())
extrap_prep{T}(::Type{Tuple{T}}, n::Val{1}) = extrap_prep(T, n)
extrap_prep{T}(::Type{Tuple{T,T}}, n::Val{1}) = extrap_prep(T, n)
extrap_prep{T}(::Type{Tuple{Tuple{T,T}}}, n::Val{1}) = extrap_prep(T, n)
function extrap_prep{S,T}(::Type{Tuple{S,T}}, n::Val{1})
    quote
        $(extrap_prep(S, n, Val{1}(), Val{:lo}()))
        $(extrap_prep(T, n, Val{1}(), Val{:hi}()))
    end
end
extrap_prep{S,T}(::Type{Tuple{Tuple{S,T}}}, n::Val{1}) = extrap_prep(Tuple{S,T}, n)

# needed for ambiguity resolution
extrap_prep{T<:Tuple}(::Type{T}, ::Val{1}) = :(throw(ArgumentError("The 1-dimensional extrap configuration $T is not supported")))

function extrap_prep{T,N}(::Type{T}, n::Val{N})
    exprs = Expr[]
    for d in 1:N
        push!(exprs, extrap_prep(T, n, Val{d}()))
    end
    return Expr(:block, exprs...)
end
function extrap_prep{N,T<:Tuple}(::Type{T}, n::Val{N})
    length(T.parameters) == N || return :(throw(ArgumentError("The $N-dimensional extrap configuration $T is not supported - must be a tuple of length $N (was length $(length(T.parameters)))")))
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
extrap_prep{T,N,d}(::Type{T}, n::Val{N}, dim::Val{d}) = extrap_prep(Tuple{T,T}, n, dim)
function extrap_prep{S,T,N,d}(::Type{Tuple{S,T}}, n::Val{N}, dim::Val{d})
    quote
        $(extrap_prep(S, n, dim, Val{:lo}()))
        $(extrap_prep(T, n, dim, Val{:hi}()))
    end
end
