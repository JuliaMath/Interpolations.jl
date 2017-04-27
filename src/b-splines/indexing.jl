using Base.Cartesian
using DualNumbers

import Base.getindex

@compat Base.IndexStyle(::Type{<:AbstractInterpolation}) = IndexCartesian()

define_indices{IT}(::Type{IT}, N, pad) = Expr(:block, Expr[define_indices_d(iextract(IT, d), d, padextract(pad, d)) for d = 1:N]...)

coefficients{IT}(::Type{IT}, N) = Expr(:block, Expr[coefficients(iextract(IT, d), N, d) for d = 1:N]...)

function gradient_coefficients{IT<:DimSpec{BSpline}}(::Type{IT}, N, dim)
    exs = Expr[d==dim ? gradient_coefficients(iextract(IT, dim), d) :
                        coefficients(iextract(IT, d), N, d) for d = 1:N]
    Expr(:block, exs...)
end
function hessian_coefficients{IT<:DimSpec{BSpline}}(::Type{IT}, N, dim1, dim2)
    exs = if dim1 == dim2
        Expr[d==dim1==dim2 ? hessian_coefficients(iextract(IT, dim), d) :
                             coefficients(iextract(IT, d), N, d) for d in 1:N]
    else
        Expr[d==dim1 || d==dim2 ? gradient_coefficients(iextract(IT, dim), d) :
                                  coefficients(iextract(IT, d), N, d) for d in 1:N]
    end
    Expr(:block, exs...)
end

index_gen{IT}(::Type{IT}, N::Integer, offsets...) = index_gen(iextract(IT, min(length(offsets)+1, N)), IT, N, offsets...)

function getindex_impl{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(itp::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}})
    meta = Expr(:meta, :inline)
    quote
        $meta
        @nexprs $N d->(x_d = xs[d])
        inds_itp = indices(itp)

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(IT, N, Pad))

        # Calculate coefficient weights based on fx
        $(coefficients(IT, N))

        # Generate the indexing expression
        @inbounds ret = $(index_gen(IT, N))
        ret
    end
end

@generated function getindex{T,N}(itp::BSplineInterpolation{T,N}, xs::Number...)
    getindex_impl(itp)
end

function gradient_impl{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(itp::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}})
    meta = Expr(:meta, :inline)
    # For each component of the gradient, alternately calculate
    # coefficients and set component
    n = count_interp_dims(IT, N)
    exs = Array{Expr, 1}(2n)
    cntr = 0
    for d = 1:N
        if count_interp_dims(iextract(IT, d), 1) > 0
            cntr += 1
            exs[2cntr-1] = gradient_coefficients(IT, N, d)
            exs[2cntr] = :(@inbounds g[$cntr] = $(index_gen(IT, N)))
        end
    end
    gradient_exprs = Expr(:block, exs...)
    quote
        $meta
        length(g) == $n || throw(ArgumentError(string("The length of the provided gradient vector (", length(g), ") did not match the number of interpolating dimensions (", n, ")")))
        @nexprs $N d->(x_d = xs[d])
        inds_itp = indices(itp)

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(IT, N, Pad))

        $gradient_exprs

        g
    end
end

# there is a Heisenbug, when Base.promote_op is inlined into getindex_return_type
# thats why we use this @noinline fence
@noinline _promote_mul(a,b) = Base.promote_op(@functorize(*), a, b)

@noinline function getindex_return_type{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}}, argtypes::Tuple)
    reduce(_promote_mul, eltype(TCoefs), argtypes)
end

function getindex_return_type{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad,I}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}}, ::Type{I})
    _promote_mul(eltype(TCoefs), I)
end

@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, xs::Number...)
    length(xs) == N || error("Can only be called with $N indexes")
    gradient_impl(itp)
end

@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, index::CartesianIndex{N})
    args = [:(index[$d]) for d = 1:N]
    :(gradient!(g, itp, $(args...)))
end

# @eval uglyness required for disambiguation with method in Base
for R in [:Real, :Any]
    @eval @generated function gradient{T,N}(itp::AbstractInterpolation{T,N}, xs::$R...)
        n = count_interp_dims(itp, N)
        Tg = promote_type(T, [x <: AbstractArray ? eltype(x) : x for x in xs]...)
        xargs = [:(xs[$d]) for d in 1:length(xs)]
        :(gradient!(Array{$Tg, 1}($n), itp, $(xargs...)))
    end
end

gradient1{T}(itp::AbstractInterpolation{T,1}, x) = gradient(itp, x)[1]

function hessian_impl{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(itp::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}})
    meta = Expr(:meta, :inline)
    # For each component of the hessian, alternately calculate
    # coefficients and set component
    n = count_interp_dims(IT, N)
    exs = Expr[]
    cntr = 0
    for d1 in 1:N, d2 in 1:N
        if count_interp_dims(iextract(IT,d1), 1) > 0 && count_interp_dims(iextract(IT,d2),1) > 0
            cntr += 1
            push!(exs, hessian_coefficients(IT, N, d1, d2))
            push!(exs, :(@inbounds H[$cntr] = $(index_gen(IT, N))))
        end
    end
    hessian_exprs = Expr(:block, exs...)

    quote
        $meta
        size(H) == ($n,$n) || throw(ArgumentError(string("The size of the provided Hessian matrix wasn't a square matrix of size ", size(H))))
        @nexprs $N d->(x_d = xs[d])
        inds_itp = indices(itp)

        $(define_indices(IT, N, Pad))

        $hessian_exprs

        H
    end
end

@generated function hessian!{T,N}(H::AbstractMatrix, itp::BSplineInterpolation{T,N}, xs::Number...)
    length(xs) == N || throw(ArgumentError("Can only be called with $N indexes"))
    hessian_impl(itp)
end

@generated function hessian!{T,N}(H::AbstractMatrix, itp::BSplineInterpolation{T,N}, index::CartesianIndex{N})
    args = [:(index[$d]) for d in 1:N]
    :(hessian!(H, itp, $(args...)))
end

@generated function hessian{T,N}(itp::AbstractInterpolation{T,N}, xs...)
    n = count_interp_dims(itp,N)
    TH = promote_type(T, [x <: AbstractArray ? eltype(x) : x for x in xs]...)
    xargs = [:(xs[$d]) for d in 1:length(xs)]
    :(hessian!(Array{TH, 2}($n,$n), itp, $(xargs...)))
end

hessian1{T}(itp::AbstractInterpolation{T,1}, x) = hessian(itp, x)[1,1]
