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

A2 = rand(Float64, nx, nx) * 100
for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
    itp_a = interpolate(A2, Tuple{BSpline(Linear), BSpline(Quadratic(BC))}, GT)
    itp_b = interpolate(A2, Tuple{BSpline(Quadratic(BC)), BSpline(Linear)}, GT)

    for i = 1:10
        x = rand()*(nx-2)+1.5
        y = rand()*(nx-2)+1.5
        xd = dual(x, 1)
        yd = dual(y, 1)
        gtmp = gradient(itp_a, x, y)
        @test_approx_eq epsilon(itp_a[xd, y]) gtmp[1]
        @test_approx_eq epsilon(itp_a[x, yd]) gtmp[2]
        gtmp = gradient(itp_b, x, y)
        @test_approx_eq epsilon(itp_b[xd, y]) gtmp[1]
        @test_approx_eq epsilon(itp_b[x, yd]) gtmp[2]
    end
end

end
