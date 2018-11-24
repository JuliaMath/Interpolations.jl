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
        if A isa AbstractRange
            show(io, A)
        else
            Base.showarg(io, A, false)
        end
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

function Base.showarg(io::IO, A::MonotonicInterpolation{T, TCoeffs, Tel, Type, K}, toplevel) where {T, TCoeffs, Tel, Type, K}
    print(io, "interpolate(")
    _showknots(io, A.knots)
    print(io, ", ")
    Base.showarg(io, A.A, false)
    print(io, ", ")
    show(io, A.it)
    if toplevel
        print(io, ") with element type ",T)
    else
        print(io, ')')
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

"""
    show_ranged(io, X, knots)

A replacement for the default array-`show` for types that may not have the canonical evaluation points.
`rngs` is the tuple of knots along each axis.
"""
function show_ranged(io::IO, X, knots)
    summary(io, X)
    isempty(X) && return
    print(io, ":")
    if !(haskey(io, :compact)) && length(axes(X, 2)) > 1
        io = IOContext(io, :compact => true)
    end
    if get(io, :limit, false) && eltype(X) === Method
        io = IOContext(io, :limit => false)
    end
    if get(io, :limit, false) && (displaysize(io))[1] - 4 <= 0
        return print(io, " …")
    else
        println(io)
    end
    io = IOContext(io, :typeinfo => eltype(X))
    Base.print_array(io, [X(x...) for x in Iterators.product(knots...)])
end

getknots(X::BSplineInterpolation)  = axes(X)
getknots(X::ScaledInterpolation)   = X.ranges
getknots(X::GriddedInterpolation)  = X.knots
getknots(X::AbstractExtrapolation) = getknots(parent(X))

Base.show(io::IO, ::MIME{Symbol("text/plain")}, X::ScaledInterpolation)   = show_ranged(io, X, getknots(X))
Base.show(io::IO, ::MIME{Symbol("text/plain")}, X::GriddedInterpolation)  = show_ranged(io, X, getknots(X))
Base.show(io::IO, ::MIME{Symbol("text/plain")}, X::AbstractExtrapolation) = show_ranged(io, X, getknots(X))
