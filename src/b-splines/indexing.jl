using Base.Cartesian

import Base.getindex

Base.linearindexing{T<:AbstractInterpolation}(::Type{T}) = Base.LinearSlow()

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
        length(g) == $n || throw(DimensionMismatch("Gradient has wrong number of components"))
        @nexprs $N d->(x_d = xs[d])

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(IT, N, Pad))

        $gradient_exprs

        g
    end
end


@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, xs::Number...)
    length(xs) == N || error("Can only be called with $N indexes")
    gradient_impl(itp)
end

@generated function gradient!{T,N}(g::AbstractVector, itp::BSplineInterpolation{T,N}, index::CartesianIndex{N})
    args = [:(index[$d]) for d = 1:N]
    :(gradient!(g, itp, $(args...)))
end

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")
