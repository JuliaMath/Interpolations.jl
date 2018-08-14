module QuadraticTests

using Interpolations, Test

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A), BSpline(Quadratic(BC())), GT()))
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

    f(x,y) = sin(x/10)*cos(y/6)
    xmax, ymax = 30,10
    A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

    # test that inner region is close to data
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp2 = @inferred(constructor(copier(A), BSpline(Quadratic(BC())), GT()))
        isfullsize = constructor == interpolate || BC==Periodic
        if isfullsize
            @test @inferred(size(itp2)) == size(A)
            @test @inferred(axes(itp2)) == axes(A)
        else
            @test @inferred(size(itp2)) == (xmax-2,ymax-2)
            @test @inferred(axes(itp2)) == (2:xmax-1,2:ymax-1)
        end

        for x in 3.1:.2:xmax-3, y in 3.1:2:ymax-3
            @test f(x,y) ≈ itp2(x,y) atol=abs(0.1 * f(x,y))
        end
    end
end

let
    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]
    itp1 = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    @test axes(itp1) == axes(A)
    for i = 1:xmax
        @test itp1(i) ≈ A[i]
    end

    f(x,y) = sin(x/10)*cos(y/6)
    xmax, ymax = 30,10
    A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]
    itp2 = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    @test axes(itp2) == axes(A)
    for j = 1:ymax, i = 1:xmax
        @test itp2(i,j) ≈ A[i,j]
    end
end

end
