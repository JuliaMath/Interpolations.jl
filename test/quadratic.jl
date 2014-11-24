
module Quadratic1DTests
println("Testing Quadratic interpolation in 1D...")
using Interpolations, Base.Test

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

A = Float64[f(x) for x in 1:xmax]

## ExtrapError

itp1 = Interpolation(A, Quadratic(ExtendInner(),OnCell()), ExtrapError())

for x in [3.1:.2:4.3]  
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

@test_throws BoundsError itp1[-3]


## ExtrapNaN

itp2 = Interpolation(A, Quadratic(ExtendInner(),OnCell()), ExtrapNaN())

for x in [3.1:.2:4.3]  
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

xlo, xhi = itp2[.9], itp2[xmax+.2]
@test isnan(xlo)
@test isnan(xhi)

# Flat/ExtrapConstant

itp3 = Interpolation(A, Quadratic(Flat(),OnCell()), ExtrapConstant())

# Check inbounds and extrap values

for x in [3.1:.2:4.3]
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

xlo, xhi = itp3[.9], itp3[xmax+.2]
@test xlo == A[1]
@test xhi == A[end]

# Check continuity
xs = [0:.1:length(A)+1]

for i in 1:length(xs)-1
    @test_approx_eq_eps itp3[xs[i]] itp3[xs[i+1]] .1
end

# Values on all data points except edges for ExtendInner

for x in 2:xmax-1
    @test_approx_eq A[x] itp1[x]
    @test_approx_eq A[x] itp2[x]
    @test_approx_eq A[x] itp3[x]
end

end
