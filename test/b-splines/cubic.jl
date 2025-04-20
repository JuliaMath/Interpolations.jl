@testset "Cubic" begin
    for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))
        isinplace = constructor == interpolate!
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
                itp1 = @inferred(constructor(copier(A), BSpline(Cubic(BC(GT())))))
                ax1 = axes(itp1)[1]
                @test Interpolations.lbounds(itp1) == (GT == OnGrid ? (first(ax1),) : (first(ax1) - 0.5,))
                @test Interpolations.ubounds(itp1) == (GT == OnGrid ? (last(ax1),) : (last(ax1) + 0.5,))
                @test_throws ArgumentError parent(itp1)
                check_axes(itp1, A, isinplace)
                check_inbounds_values(itp1, A)
                check_oob(itp1)
                can_eval_near_boundaries(itp1)

                # test that inner region is close to data
                for x in 3.1:.2:8.1
                    @test f(x) ≈ itp1(x) atol=abs(0.1 * f(x))
                end
            end

            itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC(GT())))))
            @test_throws ArgumentError parent(itp2)
            check_axes(itp2, A2, isinplace)
            check_inbounds_values(itp2, A2)
            check_oob(itp2)
            can_eval_near_boundaries(itp2)

            for x in 3.1:.2:xmax2-3, y in 3.1:2:ymax2-3
                @test f2(x,y) ≈ itp2(x,y) atol=abs(0.1 * f2(x,y))
            end
        end
        # Issue 419: Incorrect results with BSpline(Cubic(Periodic(OnCell()))))
        let BC = Periodic, GT = OnCell
            for (A, f) in ((A0, f0), (A1, f1))
                itp1 = @inferred(constructor(copier(A), BSpline(Cubic(BC(GT())))))
                @test itp1(Interpolations.lbounds(itp1)...) ≈ itp1(Interpolations.ubounds(itp1)...)
            end
        end
    end

    # Test default constructor
    itp_type = Cubic |> BSpline
    itp = interpolate(1:5, itp_type)
    @test itp(1:0.5:5) ≈ 1:0.5:5

    ix = 1:15
    k = length(ix) - 1
    f(x) = cos((x-1)*2pi/k)
    g(x) = -2pi/k * sin((x-1)*2pi/k)

    A = map(f, ix)

    for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))

        for BC in (Line, Flat, Free, Periodic), GT in (OnGrid,OnCell)

            itp = constructor(copier(A), BSpline(Cubic(BC(GT()))))
            # test that inner region is close to data
            for x in range(ix[5], stop=ix[end-4], length=100)
                @test g(x) ≈ Interpolations.gradient1(itp,x) atol=cbrt(cbrt(eps(g(x))))
            end
        end
        # Issue 419: Incorrect results with BSpline(Cubic(Periodic(OnCell()))))
        let BC = Periodic, GT = OnCell
            itp = @inferred(constructor(A[1:end-1], BSpline(Cubic(BC(GT())))))
            @test itp(Interpolations.lbounds(itp)...) ≈ itp(Interpolations.ubounds(itp)...)
        end
    end
    itp_flat_g = interpolate(A, BSpline(Cubic(Flat(OnGrid()))))
    @test Interpolations.gradient(itp_flat_g,1)[1] ≈ 0 atol=eps()
    @test Interpolations.gradient(itp_flat_g,ix[end])[1] ≈ 0 atol=eps()

    itp_flat_c = interpolate(A, BSpline(Cubic(Flat(OnCell()))))
    @test Interpolations.gradient(itp_flat_c,0.5)[1] ≈ 0 atol=eps()
    @test Interpolations.gradient(itp_flat_c,ix[end] + 0.5)[1] ≈ 0 atol=eps()

    # Can construct from coefficients
    itp = interpolate(rand(5, 9), BSpline(Cubic(Flat(OnGrid()))))
    @test itp == Interpolations.BSplineInterpolation(itp.coefs, itp.it, itp.parentaxes)
    @test_throws ArgumentError Interpolations.BSplineInterpolation(itp.coefs, itp.it, (2:3, 2:10))

    # Issue 594
    @testset "Free BC" begin
        # should be exact for polynomials of degree <= 3 (up to rounding error)
        ix = 1:15
        f = x -> x^3 - 2x^2 + 3x - 4
        A = map(f, ix)
        itp = interpolate(A, BSpline(Cubic(Free(OnGrid()))))
        @test all(x -> itp(x) ≈ f(x), 1:0.1:15)

        # also in 2d
        f2(x, y) = x^3 - 2x^2 + 3x - 4 + 2y^3 - 2y^2 + 3y - 4
        iy = 1:15
        A = map(I -> f2(I[1], I[2]), Iterators.product(ix, iy))
        itp = interpolate(A, BSpline(Cubic(Free(OnGrid()))))
        @test all(x -> itp(x...) ≈ f2(x...), Iterators.product(1:0.1:15, 1:0.1:15))
    end
end
