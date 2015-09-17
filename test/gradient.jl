module GradientTests

using Base.Test, Interpolations, DualNumbers

nx = 10
f1(x) = sin((x-3)*2pi/(nx-1) - 1)
g1(x) = 2pi/(nx-1) * cos((x-3)*2pi/(nx-1) - 1)

# Gradient of Constant should always be 0
itp1 = interpolate(Float64[f1(x) for x in 1:nx-1],
            BSpline(Constant), OnGrid)

g = Array(Float64, 1)

for x in 1:nx
    @test gradient(itp1, x)[1] == 0
    @test gradient!(g, itp1, x)[1] == 0
    @test g[1] == 0
end

# Since Linear is OnGrid in the domain, check the gradients between grid points
itp1 = interpolate(Float64[f1(x) for x in 1:nx-1],
            BSpline(Linear), OnGrid)
for x in 2.5:nx-1.5
    @test_approx_eq_eps g1(x) gradient(itp1, x)[1] abs(.1*g1(x))
    @test_approx_eq_eps g1(x) gradient!(g, itp1, x)[1] abs(.1*g1(x))
    @test_approx_eq_eps g1(x) g[1] abs(.1*g1(x))
end

for i = 1:10
    x = rand()*(nx-2)+1.5
    gtmp = gradient(itp1, x)[1]
    xd = dual(x, 1)
    @test_approx_eq epsilon(itp1[xd]) gtmp
end

# Since Quadratic is OnCell in the domain, check gradients at grid points
itp1 = interpolate(Float64[f1(x) for x in 1:nx-1],
            BSpline(Quadratic(Periodic)), OnCell)
for x in 2:nx-1
    @test_approx_eq_eps g1(x) gradient(itp1, x)[1] abs(.05*g1(x))
    @test_approx_eq_eps g1(x) gradient!(g, itp1, x)[1] abs(.05*g1(x))
    @test_approx_eq_eps g1(x) g[1] abs(.1*g1(x))
end

for i = 1:10
    x = rand()*(nx-2)+1.5
    gtmp = gradient(itp1, x)[1]
    xd = dual(x, 1)
    @test_approx_eq epsilon(itp1[xd]) gtmp
end

# For a quadratic function and quadratic interpolation, we expect an
# "exact" answer
# 1d
c = 2.3
a = 8.1
o = 1.6
qfunc = x -> a*(x-c).^2 + o
dqfunc = x -> 2*a*(x-c)
xg = Float64[1:5;]
y = qfunc(xg)

iq = interpolate(y, BSpline(Quadratic(Free)), OnCell)
x = 1.8
@test_approx_eq iq[x] qfunc(x)
@test_approx_eq gradient(iq, x)[1] dqfunc(x)

# 2d (biquadratic)
p = [(x-1.75)^2 for x = 1:7]
A = p*p'
iq = interpolate(A, BSpline(Quadratic(Free)), OnCell)
@test_approx_eq iq[4,4] (4-1.75)^4
@test_approx_eq iq[4,3] (4-1.75)^2*(3-1.75)^2
g = gradient(iq, 4, 3)
@test_approx_eq g[1] 2*(4-1.75)*(3-1.75)^2
@test_approx_eq g[2] 2*(4-1.75)^2*(3-1.75)

iq = interpolate!(copy(A), BSpline(Quadratic(InPlace)), OnCell)
@test_approx_eq iq[4,4] (4-1.75)^4
@test_approx_eq iq[4,3] (4-1.75)^2*(3-1.75)^2
g = gradient(iq, 4, 3)
@test_approx_eq_eps g[1] 2*(4-1.75)*(3-1.75)^2 0.03
@test_approx_eq_eps g[2] 2*(4-1.75)^2*(3-1.75) 0.2

# InPlaceQ is exact for an underlying quadratic
iq = interpolate!(copy(A), BSpline(Quadratic(InPlaceQ)), OnCell)
@test_approx_eq iq[4,4] (4-1.75)^4
@test_approx_eq iq[4,3] (4-1.75)^2*(3-1.75)^2
g = gradient(iq, 4, 3)
@test_approx_eq g[1] 2*(4-1.75)*(3-1.75)^2
@test_approx_eq g[2] 2*(4-1.75)^2*(3-1.75)

A2 = rand(Float64, nx, nx) * 100
for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
    itp_a = interpolate(A2, Tuple{BSpline(Linear), BSpline(Quadratic(BC))}, GT)
    itp_b = interpolate(A2, Tuple{BSpline(Quadratic(BC)), BSpline(Linear)}, GT)
    itp_c = interpolate(A2, Tuple{NoInterp, BSpline(Quadratic(BC))}, GT)
    itp_d = interpolate(A2, Tuple{BSpline(Quadratic(BC)), NoInterp}, GT)

    for i = 1:10
        x = rand()*(nx-2)+1.5
        y = rand()*(nx-2)+1.5
        xd = dual(x, 1)
        yd = dual(y, 1)
        gtmp = gradient(itp_a, x, y)
        @test length(gtmp) == 2
        @test_approx_eq epsilon(itp_a[xd, y]) gtmp[1]
        @test_approx_eq epsilon(itp_a[x, yd]) gtmp[2]
        gtmp = gradient(itp_b, x, y)
        @test length(gtmp) == 2
        @test_approx_eq epsilon(itp_b[xd, y]) gtmp[1]
        @test_approx_eq epsilon(itp_b[x, yd]) gtmp[2]
        ix, iy = round(Int, x), round(Int, y)
        gtmp = gradient(itp_c, ix, y)
        @test length(gtmp) == 1
        @test_approx_eq epsilon(itp_c[ix, yd]) gtmp[1]
        gtmp = gradient(itp_d, x, iy)
        @test length(gtmp) == 1
        @test_approx_eq epsilon(itp_d[xd, iy]) gtmp[1]
    end
end

end
