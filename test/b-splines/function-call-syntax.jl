module ExtrapFunctionCallSyntax

using Compat.Test, Interpolations, DualNumbers
using Compat: range

# Test if b-spline interpolation by function syntax yields identical results
f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10
A = Float64[f(x) for x in 1:xmax]
itpg = interpolate(A, BSpline(Linear()), OnGrid())
schemes = (Flat,Line,Free)

for T in (Cubic, Quadratic), GC in (OnGrid, OnCell)
    for etp in map(S -> @inferred(interpolate(A, BSpline(T(S())), GC())), schemes),
        x in range(1, stop=xmax, length=100)
        @test (getindex(etp, x)) == etp(x)
    end
end

for T in (Constant, Linear), GC in (OnGrid, OnCell), x in range(1, stop=xmax, length=100)
    etp = interpolate(A, BSpline(T()), GC())
    @test (getindex(etp, x)) == etp(x)
end

end
