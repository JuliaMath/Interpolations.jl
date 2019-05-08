@testset "Linear" begin
    xmax = 10
    g1(x) = sin((x-3)*2pi/(xmax-1)-1)
    f(x) = g1(x)
    A1 = Float64[f(x) for x in 1:xmax]
    fr(x) = (x^2) // 40 + 2

    ymax = 10
    g2(y) = cos(y/6)
    f(x,y) = g1(x)*g2(y)
    A2 = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

    for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))
        isinplace = constructor == interpolate!
        itp1 = @inferred(constructor(copier(A1), BSpline(Linear())))
        itp2 = @inferred(constructor(copier(A2), BSpline(Linear())))

        @test parent(itp1) === itp1.coefs
        @test Interpolations.lbounds(itp1) == (1,)
        @test Interpolations.ubounds(itp1) == (xmax,)

        for (itp, A) in ((itp1, A1), (itp2, A2))
            check_axes(itp, A, isinplace)
            check_inbounds_values(itp, A)
            check_oob(itp)
            can_eval_near_boundaries(itp)
            I = first(eachindex(itp))
            @test itp(I) == itp(Tuple(I)...)
        end

        # Just interpolation
        for x in 1:.2:xmax
            @test f(x) ≈ itp1(x) atol=abs(0.1 * f(x))
        end

        # 2D
        for x in 2.1:.2:xmax-1, y in 1.9:.2:ymax-.9
            @test ≈(f(x,y),itp2(x,y),atol=abs(0.25 * f(x,y)))
        end

        # Rational element types
        A1R = Rational{Int}[fr(x) for x in 1:10]
        itp1r = @inferred(constructor(copier(A1R), BSpline(Linear())))
        @test @inferred(size(itp1r)) == size(A1R)
        @test itp1r(23 // 10) ≈ fr(23 // 10) atol=abs(0.1 * fr(23 // 10))
        @test typeof(itp1r(23//10)) == Rational{Int}
        @test eltype(itp1r) == Rational{Int}
    end

    # Issue #183
    x = rand(3,3,3)
    itp = interpolate(x, BSpline(Linear()))
    @test itp(1.5, CartesianIndex((2, 3))) === itp(1.5, 2, 3)
    @test itp(CartesianIndex((1, 2)), 1.5) === itp(1, 2, 1.5)

    # Issue #289
    for (x, it) in ((ones(1, 3), (NoInterp(), BSpline(Linear()))),
                    (ones(3, 1), (BSpline(Linear()), NoInterp())))
        @test_throws ArgumentError interpolate(x, BSpline(Linear()))
        itp = interpolate(x, it)
        etp = extrapolate(itp, Flat())
        @test sum(isnan.(etp)) == 0
    end
end
