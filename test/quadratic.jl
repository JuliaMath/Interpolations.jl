
module Quadratic1DTests
println("Testing Quadratic interpolation in 1D...")
using Interpolations, Base.Test

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

A = Float64[f(x) for x in 1:xmax]

## ExtrapError

itp1 = Interpolation(A, Quadratic(Flat(),OnCell()), ExtrapError())

for x in [3.1:.2:4.3]  
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

@test_throws BoundsError itp1[-3]


## ExtrapNaN

itp2 = Interpolation(A, Quadratic(Flat(),OnCell()), ExtrapNaN())

for x in [3.1:.2:4.3]
    @test_approx_eq_eps f(x) itp2[x] abs(.1*f(x))
end

@test isnan(itp2[.4])
@test !isnan(itp2[.5])
@test !isnan(itp2[.6])
@test !isnan(itp2[xmax+.4])
@test !isnan(itp2[xmax+.5])
@test isnan(itp2[xmax+.6])

# Flat/ExtrapConstant

itp3 = Interpolation(A, Quadratic(Flat(),OnGrid()), ExtrapConstant())

# Check inbounds and extrap values

for x in [3.1:.2:4.3]
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

xlo, xhi = itp3[.9], itp3[xmax+.2]
@test_approx_eq xlo A[1]
@test_approx_eq xhi A[end]

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
