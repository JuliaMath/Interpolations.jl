module LinearTests

using Interpolations
using Base.Test

xmax = 10
g1(x) = sin((x-3)*2pi/(xmax-1)-1)
f(x) = g1(x)

A1 = Float64[f(x) for x in 1:xmax]

itp1c = @inferred(interpolate(A1, BSpline(Linear), OnCell))

# Just interpolation
for x in 1:.2:xmax
    @test_approx_eq_eps f(x) itp1c[x] abs(.1*f(x))
end

# Rational element types
fr(x) = (x^2) // 40 + 2
A1R = Rational{Int}[fr(x) for x in 1:10]
itp1r = @inferred(interpolate(A1R, BSpline(Linear), OnGrid))
@test_approx_eq_eps itp1r[23//10] fr(23//10) abs(.1*fr(23//10))
@test typeof(itp1r[23//10]) == Rational{Int}
@test eltype(itp1r) == Rational{Int}

# 2D
g2(y) = cos(y/6)
f(x,y) = g1(x)*g2(y)
ymax = 10
A2 = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]
itp2 = @inferred(interpolate(A2, BSpline(Linear), OnGrid))

for x in 2.1:.2:xmax-1, y in 1.9:.2:ymax-.9
    @test_approx_eq_eps f(x,y) itp2[x,y] abs(.25*f(x,y))
end

end