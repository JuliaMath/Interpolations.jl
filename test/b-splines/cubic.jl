module CubicTests

using Base.Test
using Interpolations

for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))
    f0(x) = sin((x-3)*2pi/9 - 1)
    f1(x) = 1.0 + 0.1*x + 0.01*x^2 + 0.001*x^3

    xmax = 10
    A0 = Float64[f0(x) for x in 1:xmax]
    A1 = Float64[f1(x) for x in 1:xmax]

    f2(x, y) = sin(x/10)*cos(y/6)
    xmax2, ymax2 = 30, 10
    A2 = Float64[f2(x, y) for x in 1:xmax2, y in 1:ymax2]

    for BC in (Line, Flat, Free, Periodic), GT in (OnGrid, OnCell)
        for (A, f) in ((A0, f0), (A1, f1))
            itp1 = @inferred(constructor(copier(A), BSpline(Cubic(BC())), GT()))
            @test @inferred(size(itp1)) == size(A)

            # test that inner region is close to data
            for x in 3.1:.2:8.1
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

        itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC())), GT()))
        @test @inferred(size(itp2)) == size(A2)

        for x in 3.1:.2:xmax2-3, y in 3.1:2:ymax2-3
            @test_approx_eq_eps f2(x,y) itp2[x,y] abs(.1*f2(x,y))
        end
    end
end

end

module CubicGradientTests

using Interpolations, Base.Test

ix = 1:15
f(x) = cos((x-1)*2pi/(length(ix)-1))
g(x) = -2pi/14 * sin((x-1)*2pi/(length(ix)-1))

A = f(ix)

for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))

    for BC in (Line, Flat, Free, Periodic), GT in (OnGrid,OnCell)

        itp = constructor(copier(A), BSpline(Cubic(BC())), GT())
        # test that inner region is close to data
        for x in linspace(ix[5],ix[end-4],100)
            @test_approx_eq_eps g(x) gradient(itp,x)[1] cbrt(cbrt(eps(g(x))))
        end
    end
end
itp_flat_g = interpolate(A, BSpline(Cubic(Flat())), OnGrid())
@test_approx_eq_eps gradient(itp_flat_g, 1)[1] 0 eps()
@test_approx_eq_eps gradient(itp_flat_g, ix[end])[1] 0 eps()

itp_flat_c = interpolate(A, BSpline(Cubic(Flat())), OnCell())
@test_approx_eq_eps gradient(itp_flat_c, .5)[1] 0 eps()
@test_approx_eq_eps gradient(itp_flat_c, ix[end]+.5)[1] 0 eps()

end
