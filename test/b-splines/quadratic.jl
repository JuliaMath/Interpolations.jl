@testset "Quadratic" begin
    for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
        isinplace = constructor == interpolate!
        f(x) = sin((x-3)*2pi/9 - 1)
        xmax = 10
        A = Float64[f(x) for x in 1:xmax]
        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp1 = @inferred(constructor(copier(A), BSpline(Quadratic(BC(GT())))))
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

        f(x,y) = sin(x/10)*cos(y/6)
        xmax, ymax = 30,10
        A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

        # test that inner region is close to data
        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp2 = @inferred(constructor(copier(A), BSpline(Quadratic(BC(GT())))))
            check_axes(itp2, A, isinplace)
            check_inbounds_values(itp2, A)
            check_oob(itp2)
            can_eval_near_boundaries(itp2)

            for x in 3.1:.2:xmax-3, y in 3.1:2:ymax-3
                @test f(x,y) ≈ itp2(x,y) atol=abs(0.1 * f(x,y))
            end
        end
    end

    # InPlace
    let
        f(x) = sin((x-3)*2pi/9 - 1)
        xmax = 10
        A = Float64[f(x) for x in 1:xmax]
        itp1 = interpolate!(copy(A), BSpline(Quadratic(InPlace(OnCell()))))
        @test axes(itp1) == axes(A)
        check_inbounds_values(itp1, A)

        f(x,y) = sin(x/10)*cos(y/6)
        xmax, ymax = 30,10
        A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]
        itp2 = interpolate!(copy(A), BSpline(Quadratic(InPlace(OnCell()))))
        @test axes(itp2) == axes(A)
        check_inbounds_values(itp2, A)
    end
end
