
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

xlo, xhi = itp3[.2], itp3[xmax+.7]
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

# Rational element types
A = Rational{Int}[x^2//10 for x in 1:10]
itp = Interpolation(A, Quadratic(Free(),OnCell()), ExtrapNaN())
@test itp[11//10] == (11//10)^2//10

end

module Quadratic2DTests
    println("Testing Quadratic interpolation in 2D...")
    using Interpolations, Base.Test

    f(x,y) = sin(x/10)*cos(y/6)

    xmax, ymax = 30,10

    A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

    itp1 = Interpolation(A, Quadratic(Line(), OnGrid()), ExtrapError())

    for x in 2.1:.2:xmax-1, y in 1.9:.2:ymax-.9
        @test_approx_eq_eps f(x,y) itp1[x,y] .05
    end
end
