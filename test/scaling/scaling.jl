using Interpolations
using Test, LinearAlgebra, StaticArrays

@testset "Scaling" begin
    # Model linear interpolation of y = -3 + .5x by interpolating y=x
    # and then scaling to the new x range

    itp = interpolate(1:1.0:10, BSpline(Linear()))

    sitp = @inferred(scale(itp, -3:.5:1.5))
    @test typeof(sitp) <: Interpolations.ScaledInterpolation
    @test parent(sitp) === itp

    for (x,y) in zip(-3:.05:1.5, 1:.1:10,)
        @test sitp(x) ≈ y
    end
    # trailing 1s, issue #301
    @test sitp(0.8, 1) == sitp(0.8)


    @test_throws ArgumentError scale(itp, reverse(-3:.5:1.5))

    # Model linear interpolation of y = -3 + .5x by interpolating y=x
    # and then scaling to the new x range, this time using a LinRange,
    # which should give the same values as with the StepRangeLen used previously

    itp = interpolate(1:1.0:10, BSpline(Linear()))

    sitp = @inferred(scale(itp, LinRange(-3,1.5,10)))
    @test typeof(sitp) <: Interpolations.ScaledInterpolation

    for (x,y) in zip(-3:.05:1.5, 1:.1:10,)
        @test sitp(x) ≈ y
    end

    # Verify that it works in >1D, with different types of ranges

    gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
    testfunction(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)

    xs = -5:.5:5
    ys = -4:.2:4
    zs = Float64[testfunction(x,y) for x in xs, y in ys]

    itp2 = interpolate(zs, BSpline(Quadratic(Flat(OnGrid()))))
    sitp2 = @inferred scale(itp2, xs, ys)

    for x in xs, y in ys
        @test testfunction(x,y) ≈ sitp2(x,y)
    end

    # Test gradients of scaled grids
    xs = -pi:.1:pi
    ys = map(sin, xs)
    itp = interpolate(ys, BSpline(Linear()))
    sitp = @inferred scale(itp, xs)

    for x in -pi:.1:pi
        g = @inferred(Interpolations.gradient(sitp, x))[1]
        @test cos(x) ≈ g atol=0.05
    end

    # Test Hessians of scaled grids
    xs = -pi:.1:pi
    ys = -pi:.2:pi
    zs = sin.(xs) .* sin.(ys')
    itp = interpolate(zs, BSpline(Cubic(Line(OnGrid()))))
    sitp = @inferred scale(itp, xs, ys)

    for x in xs[2:end-1], y in ys[2:end-1]
        h = @inferred(Interpolations.hessian(sitp, x, y))
        @test issymmetric(h)
        @test [-sin(x) * sin(y) cos(x) * cos(y)
                cos(x) * cos(y) -sin(x) * sin(y)] ≈ h atol=0.03
    end

    # Test Hessians of scaled grids with LinRange
    xs = LinRange(-pi, pi, 62)
    ys = LinRange(-pi, pi, 32)
    zs = sin.(xs) .* sin.(ys')
    itp = interpolate(zs, BSpline(Cubic(Line(OnGrid()))))
    sitp = @inferred scale(itp, xs, ys)

    for x in xs[2:end-1], y in ys[2:end-1]
        h = @inferred(Interpolations.hessian(sitp, x, y))
        @test issymmetric(h)
        @test [-sin(x) * sin(y) cos(x) * cos(y)
                cos(x) * cos(y) -sin(x) * sin(y)] ≈ h atol=0.03
    end

    # Verify that return types are reasonable
    @inferred(sitp2(-3.4, 1.2))
    @inferred(sitp2(-3, 1))
    @inferred(sitp2(-3.4, 1))

    sitp32 = @inferred scale(interpolate(Float32[testfunction(x,y) for x in -5:.5:5, y in -4:.2:4], BSpline(Quadratic(Flat(OnGrid())))), -5f0:.5f0:5f0, -4f0:.2f0:4f0)
    @test typeof(@inferred(sitp32(-3.4f0, 1.2f0))) == Float32

    # Iteration
    itp = interpolate(rand(3,3,3), BSpline(Quadratic(Flat(OnCell()))))
    knots = map(d->1:10:21, 1:3)
    sitp = @inferred scale(itp, knots...)

    iter = @inferred(eachvalue(sitp))

    function foo!(dest, sitp)
        i = 0
        for s in eachvalue(sitp)
            dest[i+=1] = s
        end
        dest
    end
    function bar!(dest, sitp)
        for I in CartesianIndices(size(dest))
            dest[I] = sitp(I)
        end
        dest
    end
    rfoo = Array{Float64}(undef, Interpolations.ssize(sitp))
    rbar = similar(rfoo)
    foo!(rfoo, sitp)
    bar!(rbar, sitp)
    @test rfoo ≈ rbar
end

@testset "eachvalue iteration" begin
    A = reshape(reinterpret(SVector{2,Float64}, collect(1.:12.)), 3, 2)
    sitp = scale(interpolate(A, BSpline(Cubic(Line(OnGrid())))), 1:3, 1:2)
    @test first(eachvalue(sitp)) ≈ @SVector [1.,2.]
end

@testset "axes" begin
    saxs = -50:10:50
    itp = interpolate(-5:5, BSpline(Linear()))
    sitp = @inferred(scale(itp, saxs))
    @test collect(eachvalue(sitp)) ≈ OffsetArray(-5:0.1:5, -50:50)
end
