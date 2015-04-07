using Base.Cartesian

import Base.getindex

function getindex_impl{T,N,TCoefs,IT<:BSpline,GT<:GridType}(itp::Type{BSplineInterpolation{T,N,TCoefs,IT,GT}})
    quote
        @nexprs $N d->(x_d = xs[d])

        # Calculate the indices of all coefficients that will be used
        # and define fx = x - xi in each dimension
        $(define_indices(IT, GT, N))

        # Calculate coefficient weights based on fx
        $(coefficients(IT, GT, N))

        # Generate the indexing expression
        @inbounds ret = $(index_gen(IT, N))
        ret
    end
end

# Resolve ambiguity with general array indexing,
#   getindex{T,N}(A::AbstractArray{T,N}, I::AbstractArray{T,N})
function getindex{T,N}(itp::BSplineInterpolation{T,N}, I::AbstractArray{T,N})
    error("Array indexing is not defined for interpolation objects.")
end

# Resolve ambiguity with colon indexing for 1D interpolations
#   getindex{T}(A::AbstractArray{T,1}, C::Colon)
function getindex{T}(itp::BSplineInterpolation{T,1}, c::Colon)
    error("Colon indexing is not supported for interpolation objects")
end

# Resolve ambiguity with indexing with Real indices
#   getindex{T,N}(A::AbstractArray{T,N}, x::Real...)
stagedfunction getindex{T,N,TCoefs,IT<:BSpline}(itp::BSplineInterpolation{T,N,TCoefs,IT}, xs::Real...)
    getindex_impl(itp)
end

# Linear indexing is supported only for 1D interpolations
stagedfunction getindex{T,N}(itp::BSplineInterpolation{T,N}, xs::Real)
    if N > 1
        error("Linear indexing is not supported for interpolation objects")
    end
    getindex_impl(itp)
end

stagedfunction getindex{T,N}(itp::BSplineInterpolation{T,N}, xs...)
    getindex_impl(itp)
end

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")
