
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

# Values on all data points except edges for ExtendInner

for x in 2:xmax-1
    @test_approx_eq A[x] itp1[x]
    @test_approx_eq A[x] itp2[x]
end

end
