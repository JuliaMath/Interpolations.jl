function Base.showarg(io::IO, A::BSplineInterpolation{T,N,TW,ST}, toplevel) where {T,N,TW,ST}
    print(io, "interpolate(")
    Base.showarg(io, A.coefs, false)
    print(io, ", ")
    show(io, itpflag(A))
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
    show(io, itpflag(A))
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

function Base.showarg(io::IO, A::Extrapolation{T,N,TI,IT,ET}, toplevel) where {T,N,TI,IT,ET}
    print(io, "extrapolate(")
    Base.showarg(io, A.itp, false)
    print(io, ", ")
    show(io, etpflag(A))
    print(io, ')')
    if toplevel
        print(io, " with element type ",T)
    end
end

function Base.showarg(io::IO, A::FilledExtrapolation{T,N,TI,IT}, toplevel) where {T,N,TI,IT}
    print(io, "extrapolate(")
    Base.showarg(io, A.itp, false)
    print(io, ", ", A.fillvalue, ')')
    if toplevel
        print(io, " with element type ",T)
    end
end
