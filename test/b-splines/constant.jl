@testset "Constant" begin
    # Instantiation
    N1 = 10
    A1 = rand(Float64, N1) * 100
    A2 = rand(Float64, N1, N1) * 100
    A3 = rand(Float64, N1, N1, N1) * 100

    @testset "constructor" begin
        @test Constant() == Constant{Nearest}()
        for T in (Nearest, Previous, Next)
            it = Constant{T}()
            @test it isa Constant{T}
            it = T |> Constant
            @test it isa Constant{T}
        end
    end

    @testset "Constant{Nearest}" begin
        for (constructor, copier) in ((interpolate, x -> x), (interpolate!, copy))
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
            for i in 2:N1 - 1
                @test A1[i] == itp1(i + .3) == itp1(i - .3)
            end
            # 2D
            for i in 2:N1 - 1, j in 2:N1 - 1
                @test A2[i,j] == itp2(i + .4, j - .3) == itp2(i - .4, j + .3)
            end
            # 3D
            for i in 2:N1 - 1, j in 2:N1 - 1, k in 2:N1 - 1
                @test A3[i,j,k] == itp3(i + .4, j - .3, k + .1) == itp3(i + .4, j - .3, k + .2)
            end
        end
    end

    @testset "Constant{Previous}" begin
        for (constructor, copier) in ((interpolate, x -> x), (interpolate!, copy))
            isinplace = constructor == interpolate!
            itp1 = @inferred(constructor(copier(A1), BSpline(Constant{Previous}())))
            itp2 = @inferred(constructor(copier(A2), BSpline(Constant{Previous}())))
            itp3 = @inferred(constructor(copier(A3), BSpline(Constant{Previous}())))

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
            for i in 2:N1 - 1
                @test A1[i] == itp1(i + .3) == itp1(i + .6) 
                @test A1[i - 1] == itp1(i - .3) == itp1(i - .6)
            end
            # 2D
            for i in 2:N1 - 1, j in 2:N1 - 1
                @test A2[i,j - 1] == itp2(i + .4, j - .3) == itp2(i + .4, j - .6)
            end
            # 3D
            for i in 2:N1 - 1, j in 2:N1 - 1, k in 2:N1 - 1
                @test A3[i,j - 1,k] == itp3(i + .4, j - .3, k + .1) == itp3(i + .4, j - .3, k + .2)
            end
        end
    end

    @testset "Constant{Next}" begin
        for (constructor, copier) in ((interpolate, x -> x), (interpolate!, copy))
            isinplace = constructor == interpolate!
            itp1 = @inferred(constructor(copier(A1), BSpline(Constant{Next}())))
            itp2 = @inferred(constructor(copier(A2), BSpline(Constant{Next}())))
            itp3 = @inferred(constructor(copier(A3), BSpline(Constant{Next}())))

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
            for i in 2:N1 - 1
                @test A1[i + 1] == itp1(i + .3) == itp1(i + .6) 
                @test A1[i] == itp1(i - .3) == itp1(i - .6)
            end
            # 2D
            for i in 2:N1 - 1, j in 2:N1 - 1
                @test A2[i + 1,j] == itp2(i + .4, j - .3) == itp2(i + .4, j - .6)
            end
            # 3D
            for i in 2:N1 - 1, j in 2:N1 - 1, k in 2:N1 - 1
                @test A3[i + 1,j,k + 1] == itp3(i + .4, j - .3, k + .1) == itp3(i + .4, j - .3, k + .2)
            end
        end
    end

    @testset "Constant periodic" begin
        # Constructors
        @test Constant() === Constant(Throw(OnGrid())) === Constant{Nearest,Throw{OnGrid}}()
        @test Constant() isa Constant{Nearest,Throw{OnGrid}}
        @test Constant(Periodic()) === Constant(Periodic(OnCell())) === Constant{Nearest,Periodic{OnCell}}()
        @test Constant(Periodic()) isa Constant{Nearest,Periodic{OnCell}}
        for T in (Nearest, Previous, Next)
            it = Constant{T}()
            @test it isa Constant{T, Throw{OnGrid}}
            @test "$it" == "Constant{$T}()"
            it = Constant{T}(Periodic())
            @test it isa Constant{T, Periodic{OnCell}}
            @test "$it" == "Constant{$T}(Periodic(OnCell()))"
        end

        for (constructor, copier) in ((interpolate, x -> x), (interpolate!, copy))
            isinplace = constructor == interpolate!
            itp_periodic = @inferred(constructor(copier(A1), BSpline(Constant(Periodic()))))
            itp_previous = @inferred(constructor(copier(A1), BSpline(Constant{Previous}(Periodic()))))
            itp_next = @inferred(constructor(copier(A1), BSpline(Constant{Next}(Periodic()))))

            for itp in (itp_periodic, itp_previous, itp_next)
                @test parent(itp) === itp.coefs
                @test all(Interpolations.lbounds(itp) .≈ (0.5,))
                @test all(Interpolations.ubounds(itp) .≈ (N1 + 0.5,))

                check_axes(itp, A1, isinplace)
                check_inbounds_values(itp, A1)
                check_oob(itp)
                can_eval_near_boundaries(itp)
            end

            # Evaluation between data points (tests constancy)
            for i in 2:N1-1
                @test A1[i] == itp_periodic(i + .3) == itp_periodic(i - .3)
                @test A1[i] == itp_previous(i + .3) == itp_previous(i + .6)
                @test A1[i - 1] == itp_previous(i - .3) == itp_previous(i - .6)
                @test A1[i + 1] == itp_next(i + .3) == itp_next(i + .6)
                @test A1[i] == itp_next(i - .3) == itp_next(i - .6)
            end

            # Evaluation between data points in [0.5, 1.5], [N1 - 0.5, N1 + 0.5].
            @test A1[1] == itp_periodic(1 - .3) == itp_periodic(1 + .3)
            @test A1[N1] == itp_periodic(N1 - .3) == itp_periodic(N1 + .3)
            @test A1[1] == itp_previous(1 + .3)
            @test A1[N1] == itp_previous(N1 + .3) == itp_previous(1 - .3)
            @test A1[1] == itp_next(1 - .3) == itp_next(N1 + .3)
            @test A1[N1] == itp_next(N1 - .3)
        end
    end
end
