module Quadratic1DTests

using Interpolations, Base.Test

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

A = Float64[f(x) for x in 1:xmax]

for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
    itp1 = interpolate(A, BSpline(Quadratic(BC)), GT)

    # test that inner region is close to data
    for x in 3.1:.2:8.3
        @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
    end

    # test that we can evaluate close to, and at, boundaries
    if GT == OnGrid
        itp1[1.]
        itp1[1.0]
        itp1[1.2]
        itp1[9.8]
        itp1[10.]
        itp1[10]
    else
        itp1[0.5]
        itp1[0.6]
        itp1[10.4]
        itp1[10.5]
    end
end


f(x,y) = sin(x/10)*cos(y/6)
xmax, ymax = 30,10
A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

# test that inner region is close to data
for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
    itp2 = interpolate(A, BSpline(Quadratic(BC)), GT)

    for x in 3.1:.2:xmax-3, y in 3.1:2:ymax-3
        @test_approx_eq_eps f(x,y) itp2[x,y] abs(.1*f(x,y))
    end
end

end
