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
        itp1c = @inferred(constructor(copier(A1), BSpline(Linear()), OnCell()))
        itp1g = @inferred(constructor(copier(A1), BSpline(Linear()), OnGrid()))
        itp2c = @inferred(constructor(copier(A2), BSpline(Linear()), OnCell()))
        itp2g = @inferred(constructor(copier(A2), BSpline(Linear()), OnGrid()))

        @test parent(itp1c) === itp1c.coefs
        @test Interpolations.lbounds(itp1c) == (0.5,)
        @test Interpolations.lbounds(itp1g) == (1,)
        @test Interpolations.ubounds(itp1c) == (xmax + 0.5,)
        @test Interpolations.ubounds(itp1g) == (xmax,)

        for (itp, A) in ((itp1c, A1), (itp1g, A1), (itp2c, A2), (itp2g, A2))
            check_axes(itp, A, isinplace)
            check_inbounds_values(itp, A)
            check_oob(itp)
            can_eval_near_boundaries(itp)
            I = first(eachindex(itp))
            @test itp(I) == itp(Tuple(I)...)
        end

        # Just interpolation
        for x in 1:.2:xmax
            @test f(x) ≈ itp1c(x) atol=abs(0.1 * f(x))
        end

        # 2D
        for itp2 in (itp2c, itp2g)
            for x in 2.1:.2:xmax-1, y in 1.9:.2:ymax-.9
                @test ≈(f(x,y),itp2(x,y),atol=abs(0.25 * f(x,y)))
            end
        end

        # Rational element types
        A1R = Rational{Int}[fr(x) for x in 1:10]
        itp1r = @inferred(constructor(copier(A1R), BSpline(Linear()), OnGrid()))
        @test @inferred(size(itp1r)) == size(A1R)
        @test itp1r(23 // 10) ≈ fr(23 // 10) atol=abs(0.1 * fr(23 // 10))
        @test typeof(itp1r(23//10)) == Rational{Int}
        @test eltype(itp1r) == Rational{Int}
    end

    # Issue #183
    x = rand(3,3,3)
    itp = interpolate(x, BSpline(Linear()), OnGrid())
    @test itp(1.5, CartesianIndex((2, 3))) === itp(1.5, 2, 3)
    @test itp(CartesianIndex((1, 2)), 1.5) === itp(1, 2, 1.5)
end
