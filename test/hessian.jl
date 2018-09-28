using Test, Interpolations, LinearAlgebra, ForwardDiff

@testset "Hessians" begin
    nx = 5
    k = 2pi/(nx-1)
    f1(x) = sin(k*(x-3) - 1)
    A1 = Float64[f1(x) for x in 1:nx]
    h1 = Array{Float64}(undef, 1, 1)
    A2 = rand(Float64, nx, nx) * 100
    h2 = Array{Float64}(undef, 2, 2)

    for (A, h) in ((A1, h1), (A2, h2))
        for itp in (interpolate(A, BSpline(Constant())),
                    interpolate(A, BSpline(Linear())))
            if ndims(A) == 1
                # Hessian of Constant and Linear should always be 0 in 1d
                for x in InterpolationTestUtils.thirds(axes(A))
                    @test all(iszero, Interpolations.hessian(itp, x...))
                    @test all(iszero, Interpolations.hessian!(h, itp, x...))
                end
            else
                for x in InterpolationTestUtils.thirds(axes(A))
                    check_hessian(itp, h)
                end
            end
        end

        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp = interpolate(A, BSpline(Quadratic(BC(GT()))))
            check_hessian(itp, h)
            I = first(eachindex(itp))
            @test Interpolations.hessian(itp, I) == Interpolations.hessian(itp, Tuple(I)...)
        end

        for BC in (Line, Flat, Free, Periodic), GT in (OnGrid, OnCell)
            itp = interpolate(A, BSpline(Cubic(BC(GT()))))
            check_hessian(itp, h)
        end
    end

    itp = interpolate(A2, (BSpline(Quadratic(Flat(OnCell()))), NoInterp()))
    v = A2[:, 2]
    itpcol = interpolate(v, BSpline(Quadratic(Flat(OnCell()))))
    @test Interpolations.hessian(itp, 3.2, 2) == Interpolations.hessian(itpcol, 3.2)


    @testset "Monotonic" begin
        x = [0.0, 0.2, 0.5, 0.6, 0.9, 1.0]
        ys = [[-3.0, 0.0, 5.0, 10.0, 18.0, 22.0],
              [10.0, 0.0, -5.0, 10.0, -8.0, -2.0]]
        grid = 0.0:0.1:1.0

        itypes = [LinearMonotonicInterpolation(),
            FiniteDifferenceMonotonicInterpolation(),
            CardinalMonotonicInterpolation(0.0),
            CardinalMonotonicInterpolation(0.5),
            CardinalMonotonicInterpolation(1.0),
            FritschCarlsonMonotonicInterpolation(),
            FritschButlandMonotonicInterpolation(),
            SteffenMonotonicInterpolation()]

        for y in ys
            for it in itypes
                itp = interpolate(x, y, it)
                for t in grid
                    hessval = ForwardDiff.hessian(u -> itp(u[1]), [t])[1, 1]
                    @test Interpolations.hessian1(itp, t) ≈ hessval atol = 1.e-12
                    @test Interpolations.hessian(itp, t)[1] ≈ hessval atol = 1.e-12
                end
            end
        end
    end
end
