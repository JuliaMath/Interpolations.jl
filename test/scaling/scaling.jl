module ScalingTests

using Interpolations
using Base.Test

# Model linear interpolation of y = -3 + .5x by interpolating y=x
# and then scaling to the new x range

itp = interpolate(1:1.0:10, BSpline(Linear), OnGrid)

sitp = scale(itp, -3:.5:1.5)

for (x,y) in zip(-3:.05:1.5, 1:.1:10)
    @test_approx_eq sitp[x] y
end

# Verify that it works in >1D, with different types of ranges

gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
testfunction(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)

xs = -5:.5:5
ys = -4:.2:4
zs = Float64[testfunction(x,y) for x in xs, y in ys]

itp2 = interpolate(zs, BSpline(Quadratic(Flat)), OnGrid)
sitp2 = scale(itp2, xs, ys)

for x in xs, y in ys
    @test_approx_eq testfunction(x,y) sitp2[x,y]
end

# Test gradients of scaled grids
xs = -pi:.1:pi
ys = sin(xs)
itp = interpolate(ys, BSpline(Linear), OnGrid)
sitp = scale(itp, xs)

for x in -pi:.1:pi
    g = @inferred(gradient(sitp, x))[1]
    @test_approx_eq_eps cos(x) g .05
end

# Verify that return types are reasonable
@inferred(getindex(sitp2, -3.4, 1.2))
@inferred(getindex(sitp2, -3, 1))
@inferred(getindex(sitp2, -3.4, 1))

sitp32 = scale(interpolate(Float32[testfunction(x,y) for x in -5:.5:5, y in -4:.2:4], BSpline(Quadratic(Flat)), OnGrid), -5f0:.5f0:5f0, -4f0:.2f0:4f0)
@test typeof(@inferred(getindex(sitp32, -3.4f0, 1.2f0))) == Float32

end
