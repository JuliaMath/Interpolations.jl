module QuadraticTests

using Interpolations, Base.Test

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A), BSpline(Quadratic(BC())), GT()))
        @test @inferred(size(itp1)) == size(A)
        @test_throws ArgumentError parent(itp1)

        # test that inner region is close to data
        for x in 3.1:.2:8.1
            @test ≈(f(x),itp1[x],atol=abs(0.1 * f(x)))
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
        itp2 = @inferred(constructor(copier(A), BSpline(Quadratic(BC())), GT()))
        @test @inferred(size(itp2)) == size(A)

        for x in 3.1:.2:xmax-3, y in 3.1:2:ymax-3
            @test ≈(f(x,y),itp2[x,y],atol=abs(0.1 * f(x,y)))
        end
    end
end

let
    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]
    itp1 = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    for i = 1:xmax
        @test itp1[i] ≈ A[i]
    end

    f(x,y) = sin(x/10)*cos(y/6)
    xmax, ymax = 30,10
    A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]
    itp2 = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    for j = 1:ymax, i = 1:xmax
        @test itp2[i,j] ≈ A[i,j]
    end
end

end
