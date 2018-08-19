@testset "Constant" begin
    # Instantiation
    N1 = 10
    A1 = rand(Float64, N1) * 100
    A2 = rand(Float64, N1, N1) * 100
    A3 = rand(Float64, N1, N1, N1) * 100

    for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
        isinplace = constructor == interpolate!
        itp1c = @inferred(constructor(copier(A1), BSpline(Constant()), OnCell()))
        itp1g = @inferred(constructor(copier(A1), BSpline(Constant()), OnGrid()))
        itp2c = @inferred(constructor(copier(A2), BSpline(Constant()), OnCell()))
        itp2g = @inferred(constructor(copier(A2), BSpline(Constant()), OnGrid()))
        itp3c = @inferred(constructor(copier(A3), BSpline(Constant()), OnCell()))
        itp3g = @inferred(constructor(copier(A3), BSpline(Constant()), OnGrid()))

        @test parent(itp1c) === itp1c.coefs
        @test Interpolations.lbounds(itp1c) == (0.5,)
        @test Interpolations.lbounds(itp1g) == (1,)
        @test Interpolations.ubounds(itp1c) == (N1 + 0.5,)
        @test Interpolations.ubounds(itp1g) == (N1,)

        # Evaluation on provided data points
        for (itp, A) in ((itp1c, A1), (itp1g, A1), (itp2c, A2), (itp2g, A2), (itp3c, A3), (itp3g, A3))
            check_axes(itp, A, isinplace)
            check_inbounds_values(itp, A)
            check_oob(itp)
            can_eval_near_boundaries(itp)
        end

        # Evaluation between data points (tests constancy)
        for i in 2:N1-1
            @test A1[i] == itp1c(i+.3) == itp1g(i+.3) == itp1c(i-.3) == itp1g(i-.3)
        end
        # 2D
        for i in 2:N1-1, j in 2:N1-1
            @test A2[i,j] == itp2c(i+.4,j-.3) == itp2g(i+.4,j-.3)
        end
        # 3D
        for i in 2:N1-1, j in 2:N1-1, k in 2:N1-1
            @test A3[i,j,k] == itp3c(i+.4,j-.3,k+.1) == itp3g(i+.4,j-.3,k+.2)
        end

        # Edge behavior
        @test A1[1] == itp1c(.7)
        @test A1[N1] == itp1c(N1+.3)
    end
end
