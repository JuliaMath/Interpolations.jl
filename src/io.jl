using ShowItLikeYouBuildIt

Base.summary(A::AbstractInterpolation) = summary_build(A)

function ShowItLikeYouBuildIt.showarg{T,N,TW,ST,GT}(io::IO, A::BSplineInterpolation{T,N,TW,ST,GT})
    print(io, "interpolate(")
    showarg(io, A.coefs)
    print(io, ", ")
    _showtypeparam(io, ST)
    print(io, ", ")
    _showtypeparam(io, GT)
    print(io, ')')
end

function ShowItLikeYouBuildIt.showarg{T,N,TC,ST,K}(io::IO, A::GriddedInterpolation{T,N,TC,ST,K})
    print(io, "interpolate(")
    _showknots(io, A.knots)
    print(io, ", ")
    showarg(io, A.coefs)
    print(io, ", ")
    _showtypeparam(io, ST)
    print(io, ')')
end

_showknots(io, A) = showarg(io, A)
function _showknots{N}(io, tup::NTuple{N,Any})
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

function ShowItLikeYouBuildIt.showarg{T,N,TI,IT,GT,ET}(io::IO, A::Extrapolation{T,N,TI,IT,GT,ET})
    print(io, "extrapolate(")
    showarg(io, A.itp)
    print(io, ", ")
    _showtypeparam(io, ET)
    print(io, ')')
end

function ShowItLikeYouBuildIt.showarg{T,N,TI,IT,GT}(io::IO, A::FilledExtrapolation{T,N,TI,IT,GT})
    print(io, "extrapolate(")
    showarg(io, A.itp)
    print(io, ", ", A.fillvalue, ')')
end

_showtypeparam{T}(io, ::Type{T}) =
    print(io, T.name.name, "()")
_showtypeparam{T}(io, ::Type{Quadratic{T}}) =
    print(io, "Quadratic(", T.name.name, "())")
_showtypeparam{T}(io, ::Type{Cubic{T}}) =
    print(io, "Cubic(",     T.name.name, "())")

function _showtypeparam{T}(io, ::Type{BSpline{T}})
    print(io, "BSpline(")
    _showtypeparam(io, T)
    print(io, ')')
end

function _showtypeparam{T}(io, ::Type{Gridded{T}})
    print(io, "Gridded(")
    _showtypeparam(io, T)
    print(io, ')')
end

function _showtypeparam{TTup<:Tuple}(io, types::Type{TTup})
    print(io, '(')
    N = length(types.types)
    for (i, T) in enumerate(types.types)
        _showtypeparam(io, T)
        i < N && print(io, ", ")
    end
    print(io, ')')
end
