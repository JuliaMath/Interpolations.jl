using Test, Interpolations, LinearAlgebra

@testset "Hessians" begin
    nx = 5
    k = 2pi/(nx-1)
    f1(x) = sin(k*(x-3) - 1)
    A1 = Float64[f1(x) for x in 1:nx]
    h1 = Array{Float64}(undef, 1, 1)
    A2 = rand(Float64, nx, nx) * 100
    h2 = Array{Float64}(undef, 2, 2)

    for (A, h) in ((A1, h1), (A2, h2))
        for GT in (OnGrid, OnCell)
            for itp in (interpolate(A, BSpline(Constant()), GT()),
                        interpolate(A, BSpline(Linear()), GT()))
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
        end

        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp = interpolate(A, BSpline(Quadratic(BC())), GT())
            check_hessian(itp, h)
            I = first(eachindex(itp))
            @test Interpolations.hessian(itp, I) == Interpolations.hessian(itp, Tuple(I)...)
        end

        for BC in (Line, Flat, Free, Periodic), GT in (OnGrid, OnCell)
            itp = interpolate(A, BSpline(Cubic(BC())), GT())
            check_hessian(itp, h)
        end
    end

    # TODO: mixed interpolation (see gradient.jl)
end
