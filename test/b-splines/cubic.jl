module CubicTests

using Test
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
            isfullsize = constructor == interpolate || BC==Periodic
            if isfullsize
                @test @inferred(size(itp1)) == size(A)
                @test @inferred(axes(itp1)) == axes(A)
            else
                @test @inferred(size(itp1)) == (xmax-2,)
                @test @inferred(axes(itp1)) == (2:xmax-1,)
            end
            @test_throws ArgumentError parent(itp1)
            # Test that within the axes, we reconstruct exactly
            for i in eachindex(itp1)
                @test A[i] ≈ itp1[i] == itp1[Tuple(i)...] == itp1(i) ≈ itp1(float.(Tuple(i))...)
            end

            # test that inner region is close to data
            for x in 3.1:.2:8.1
                @test f(x) ≈ itp1(x) atol=abs(0.1 * f(x))
            end

            # test that we can evaluate close to, and at, boundaries
            if GT == OnGrid
                isfullsize ? @test(itp1[1] ≈ A[1]) : @test_throws BoundsError itp1[1]
                isfullsize ? @test(itp1(1.0) ≈ A[1]) : @test_throws BoundsError itp1(1.0)
                isfullsize ? itp1(1.2) : @test_throws BoundsError itp1(1.2)
                itp1(1.6)
                itp1(9.4)
                isfullsize ? itp1(9.8) : @test_throws BoundsError itp1(9.8)
                isfullsize ? @test(itp1(10.0) ≈ A[10]) : @test_throws BoundsError itp1(10.0)
                isfullsize ? @test(itp1[10] ≈ A[10]) : @test_throws BoundsError itp1[10]
            else
                isfullsize ? itp1(0.5) : @test_throws BoundsError itp1(0.5)
                isfullsize ? itp1(0.6) : @test_throws BoundsError itp1(0.6)
                itp1(1.6)
                itp1(9.4)
                isfullsize ? itp1(10.4) : @test_throws BoundsError itp1(10.4)
                isfullsize ? itp1(10.5) : @test_throws BoundsError itp1(10.5)
            end
        end

        isfullsize = constructor == interpolate || BC==Periodic
        itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC())), GT()))
        if isfullsize
            @test @inferred(size(itp2)) == size(A2)
            @test @inferred(axes(itp2)) == axes(A2)
        else
            @test @inferred(size(itp2)) == (xmax2-2,ymax2-2)
            @test @inferred(axes(itp2)) == (2:xmax2-1,2:ymax2-1)
        end

        for x in 3.1:.2:xmax2-3, y in 3.1:2:ymax2-3
            @test f2(x,y) ≈ itp2(x,y) atol=abs(0.1 * f2(x,y))
        end
    end
end

end

module CubicGradientTests

using Interpolations, Test, LinearAlgebra

ix = 1:15
k = length(ix) - 1
f(x) = cos((x-1)*2pi/k)
g(x) = -2pi/k * sin((x-1)*2pi/k)

A = map(f, ix)

for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))

    for BC in (Line, Flat, Free, Periodic), GT in (OnGrid,OnCell)

        itp = constructor(copier(A), BSpline(Cubic(BC())), GT())
        # test that inner region is close to data
        for x in range(ix[5], stop=ix[end-4], length=100)
            @test g(x) ≈ Interpolations.gradient1(itp,x) atol=cbrt(cbrt(eps(g(x))))
        end
    end
end
itp_flat_g = interpolate(A, BSpline(Cubic(Flat())), OnGrid())
@test Interpolations.gradient(itp_flat_g,1)[1] ≈ 0 atol=eps()
@test Interpolations.gradient(itp_flat_g,ix[end])[1] ≈ 0 atol=eps()

itp_flat_c = interpolate(A, BSpline(Cubic(Flat())), OnCell())
@test_broken Interpolations.gradient(itp_flat_c,0.5)[1] ≈ 0 atol=eps()
@test_broken Interpolations.gradient(itp_flat_c,ix[end] + 0.5)[1] ≈ 0 atol=eps()

end
