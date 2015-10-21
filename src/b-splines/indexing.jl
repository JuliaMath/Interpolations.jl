using Base.Cartesian

import Base.getindex

Base.linearindexing{T<:AbstractInterpolation}(::Type{T}) = Base.LinearSlow()

define_indices{IT}(::Type{IT}, N, pad) = Expr(:block, Expr[define_indices_d(iextract(IT, d), d, padextract(pad, d)) for d = 1:N]...)

coefficients{IT}(::Type{IT}, N) = Expr(:block, Expr[coefficients(iextract(IT, d), N, d) for d = 1:N]...)

function gradient_coefficients{IT<:DimSpec{BSpline}}(::Type{IT}, N, dim)
    exs = Expr[d==dim ? gradient_coefficients(iextract(IT, dim), d) :
                        coefficients(iextract(IT, d), N, d) for d = 1:N]
    Expr(:block, exs...)
end

index_gen{IT}(::Type{IT}, N::Integer, offsets...) = index_gen(iextract(IT, min(length(offsets)+1, N)), IT, N, offsets...)

function getindex_impl{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(itp::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}})
    meta = Expr(:meta, :inline)
    quote
        $meta
        @nexprs $N d->(x_d = xs[d])

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
    exs = Array(Expr, 2n)
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

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(IT, N, Pad))

        $gradient_exprs

        g
    end
end

function getindex_return_type{T,N,TCoefs,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Pad}(::Type{BSplineInterpolation{T,N,TCoefs,IT,GT,Pad}}, argtypes)
    Tret = TCoefs
    for a in argtypes
        Tret = Base.promote_op(Base.MulFun, Tret, a)
    end
    Tret
end


@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, xs::Number...)
    length(xs) == N || error("Can only be called with $N indexes")
    gradient_impl(itp)
end

@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, index::CartesianIndex{N})
    args = [:(index[$d]) for d = 1:N]
    :(gradient!(g, itp, $(args...)))
end

offsetsym(off, d) = off == -1 ? symbol("ixm_", d) :
                    off ==  0 ? symbol("ix_", d) :
                    off ==  1 ? symbol("ixp_", d) :
                    off ==  2 ? symbol("ixpp_", d) : error("offset $off not recognized")
