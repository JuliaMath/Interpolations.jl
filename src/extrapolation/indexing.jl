function extrap_impl{T<:AbstractExtrapolation}(exp::Type{T}, xs...)
    quote
        $(extrap_prep(exp, xs...))
        exp.itp[xs...]
    end
end


# Resolve ambiguity with general array indexing,
#   getindex{T,N}(A::AbstractArray{T,N}, I::AbstractArray{T,N})
function getindex{T,N}(exp::AbstractExtrapolation{T,N}, I::AbstractArray{T,N})
    error("Array indexing is not defined for interpolation objects.")
end

# Resolve ambiguity with colon indexing for 1D interpolations
#   getindex{T}(A::AbstractArray{T,1}, C::Colon)
function getindex{T}(exp::AbstractExtrapolation{T,1}, c::Colon)
    error("Colon indexing is not supported for interpolation objects")
end

# Resolve ambiguity with indexing with Real indices
#   getindex{T,N}(A::AbstractArray{T,N}, x::Real...)
@generated function getindex{T,N,ITP,GT}(exp::AbstractExtrapolation{T,N,ITP,GT}, xs::Real...)
    :($extrap_impl(exp, xs...))
end

# Linear indexing is supported only for 1D interpolations
@generated function getindex{T,N}(exp::AbstractExtrapolation{T,N}, xs::Real)
    if N > 1
        error("Linear indexing is not supported for interpolation objects")
    end
    :($(extrap_impl(exp, xs)))
end

@generated function getindex{T,N,ITP,GT}(exp::AbstractExtrapolation{T,N,ITP,GT}, xs...)
    :($(extrap_impl(exp, xs...)))
end
