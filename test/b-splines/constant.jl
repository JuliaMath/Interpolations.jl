@testset "Constant" begin
    # Instantiation
    N1 = 10
    A1 = rand(Float64, N1) * 100
    A2 = rand(Float64, N1, N1) * 100
    A3 = rand(Float64, N1, N1, N1) * 100

    for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
        isinplace = constructor == interpolate!
        itp1 = @inferred(constructor(copier(A1), BSpline(Constant())))
        itp2 = @inferred(constructor(copier(A2), BSpline(Constant())))
        itp3 = @inferred(constructor(copier(A3), BSpline(Constant())))

        @test parent(itp1) === itp1.coefs
        @test Interpolations.lbounds(itp1) == (1,)
        @test Interpolations.ubounds(itp1) == (N1,)

        # Evaluation on provided data points
        for (itp, A) in ((itp1, A1), (itp2, A2), (itp3, A3))
            check_axes(itp, A, isinplace)
            check_inbounds_values(itp, A)
            check_oob(itp)
            can_eval_near_boundaries(itp)
        end

        # Evaluation between data points (tests constancy)
        for i in 2:N1-1
            @test A1[i] == itp1(i+.3) == itp1(i+.3) == itp1(i-.3) == itp1(i-.3)
        end
        # 2D
        for i in 2:N1-1, j in 2:N1-1
            @test A2[i,j] == itp2(i+.4,j-.3) == itp2(i+.4,j-.3)
        end
        # 3D
        for i in 2:N1-1, j in 2:N1-1, k in 2:N1-1
            @test A3[i,j,k] == itp3(i+.4,j-.3,k+.1) == itp3(i+.4,j-.3,k+.2)
        end
    end
end
