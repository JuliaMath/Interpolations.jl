# after this, functionality was incorporated into Base
if VERSION < v"0.7.0-DEV.1790"
using ShowItLikeYouBuildIt

Base.summary(A::AbstractInterpolation) = summary_build(A)

function ShowItLikeYouBuildIt.showarg(io::IO, A::BSplineInterpolation{T,N,TW,ST,GT}) where {T,N,TW,ST,GT}
    print(io, "interpolate(")
    showarg(io, A.coefs)
    print(io, ", ")
    _showtypeparam(io, ST)
    print(io, ", ")
    _showtypeparam(io, GT)
    print(io, ')')
end

function ShowItLikeYouBuildIt.showarg(io::IO, A::GriddedInterpolation{T,N,TC,ST,K}) where {T,N,TC,ST,K}
    print(io, "interpolate(")
    _showknots(io, A.knots)
    print(io, ", ")
    showarg(io, A.coefs)
    print(io, ", ")
    _showtypeparam(io, ST)
    print(io, ')')
end

_showknots(io, A) = showarg(io, A)
function _showknots(io, tup::NTuple{N,Any}) where N
    print(io, '(')
    for (i, A) in enumerate(tup)
        showarg(io, A)
        i < N && print(io, ',')
    end
    N == 1 && print(io, ',')
    print(io, ')')
end

function ShowItLikeYouBuildIt.showarg(io::IO, A::ScaledInterpolation)
    print(io, "scale(")
    showarg(io, A.itp)
    print(io, ", ", A.ranges, ')')
end

function ShowItLikeYouBuildIt.showarg(io::IO, A::Extrapolation{T,N,TI,IT,GT,ET}) where {T,N,TI,IT,GT,ET}
    print(io, "extrapolate(")
    showarg(io, A.itp)
    print(io, ", ")
    _showtypeparam(io, ET)
    print(io, ')')
end

function ShowItLikeYouBuildIt.showarg(io::IO, A::FilledExtrapolation{T,N,TI,IT,GT}) where {T,N,TI,IT,GT}
    print(io, "extrapolate(")
    showarg(io, A.itp)
    print(io, ", ", A.fillvalue, ')')
end

else

function Base.showarg(io::IO, A::BSplineInterpolation{T,N,TW,ST,GT}, toplevel) where {T,N,TW,ST,GT}
    print(io, "interpolate(")
    Base.showarg(io, A.coefs, false)
    print(io, ", ")
    _showtypeparam(io, ST)
    print(io, ", ")
    _showtypeparam(io, GT)
    if toplevel
        print(io, ") with element type ",T)
    else
        print(io, ')')
    end
end

function Base.showarg(io::IO, A::GriddedInterpolation{T,N,TC,ST,K}, toplevel) where {T,N,TC,ST,K}
    print(io, "interpolate(")
    _showknots(io, A.knots)
    print(io, ", ")
    Base.showarg(io, A.coefs, false)
    print(io, ", ")
    _showtypeparam(io, ST)
    if toplevel
        print(io, ") with element type ",T)
    else
        print(io, ')')
    end
end

_showknots(io, A) = Base.showarg(io, A, false)
function _showknots(io, tup::NTuple{N,Any}) where N
    print(io, '(')
    for (i, A) in enumerate(tup)
        Base.showarg(io, A, false)
        i < N && print(io, ',')
    end
    N == 1 && print(io, ',')
    print(io, ')')
end

function Base.showarg(io::IO, A::ScaledInterpolation{T}, toplevel) where {T}
    print(io, "scale(")
    Base.showarg(io, A.itp, false)
    print(io, ", ", A.ranges, ')')
    if toplevel
        print(io, " with element type ",T)
    end
end

function Base.showarg(io::IO, A::Extrapolation{T,N,TI,IT,GT,ET}, toplevel) where {T,N,TI,IT,GT,ET}
    print(io, "extrapolate(")
    Base.showarg(io, A.itp, false)
    print(io, ", ")
    _showtypeparam(io, ET)
    print(io, ')')
    if toplevel
        print(io, " with element type ",T)
    end
end

function Base.showarg(io::IO, A::FilledExtrapolation{T,N,TI,IT,GT}, toplevel) where {T,N,TI,IT,GT}
    print(io, "extrapolate(")
    Base.showarg(io, A.itp, false)
    print(io, ", ", A.fillvalue, ')')
    if toplevel
        print(io, " with element type ",T)
    end
end

end

_showtypeparam(io, ::Type{T}) where {T} =
    print(io, T.name.name, "()")
_showtypeparam(io, ::Type{Quadratic{T}}) where {T} =
    print(io, "Quadratic(", T.name.name, "())")
_showtypeparam(io, ::Type{Cubic{T}}) where {T} =
    print(io, "Cubic(",     T.name.name, "())")

function _showtypeparam(io, ::Type{BSpline{T}}) where T
    print(io, "BSpline(")
    _showtypeparam(io, T)
    print(io, ')')
end

function _showtypeparam(io, ::Type{Gridded{T}}) where T
    print(io, "Gridded(")
    _showtypeparam(io, T)
    print(io, ')')
end

function _showtypeparam(io, types::Type{TTup}) where TTup<:Tuple
    print(io, '(')
    N = length(types.types)
    for (i, T) in enumerate(types.types)
        _showtypeparam(io, T)
        i < N && print(io, ", ")
    end
    print(io, ')')
end
