using Test, Interpolations, DualNumbers, Unitful

@testset "ExtrapTypeStability" begin

    # Test type-stability of 1-dimensional extrapolation
    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]
    itpg = interpolate(A, BSpline(Linear()))

    schemes = (
        Flat,
        Line,
        Reflect,
        Periodic
    )

    for etp in map(E -> @inferred(extrapolate(itpg, E())), schemes),
        x in (
            # In-bounds evaluation
            3.4, 3, dual(3.1),
            # Out-of-bounds evaluation
            -3.4, -3, dual(-3,1),
            13.4, 13, dual(13,1)
        )
        @inferred(etp(x))
    end

    # Test type-stability of 2-dimensional extrapolation with homogeneous scheme
    g(y) = (y/100)^3
    ymax = 4
    A = Float64[f(x)*g(y) for x in 1:xmax, y in 1:ymax]
    itp2 = interpolate(A, BSpline(Linear()))

    for (etp2,E) in map(E -> (extrapolate(itp2, E()), E), schemes),
        x in (
            # In-bounds evaluation
            3.4, 3, dual(3.1),
            # Out-of-bounds evaluation
            -3.4, -3, dual(-3,1),
            13.4, 13, dual(13,1)
        ),
        y in (
            # In-bounds evaluation
            2.1, 2, dual(2.3, 1),
            # Out-of-bounds evaluation
            -2.1, -2, dual(-2.3, 1),
            12.1, 12, dual(12.1, 1)
        )
        @inferred(etp2(x, y))
    end

    A = [1 2; 3 4]
    Af = Float64.(A)
    for B in (A, Af)
        itpg2 = interpolate(B, BSpline(Linear()))
        etp = extrapolate(itpg2, NaN)
        @test typeof(@inferred(etp(dual(1.5,1), dual(1.5,1)))) ==
              typeof(@inferred(etp(dual(6.5,1), dual(3.5,1))))
    end

    # issue #258
    t = [2017, 2018]
    u = u"V"
    y = [1.f0, 2.f0] .* u
    how = Interpolations.Gridded(Interpolations.Linear())
    itp = interpolate((t,), y, Gridded(Linear()))
    etp = extrapolate(itp, 1.f0 * u)
    @test @inferred(itp(2018)) === 2.0f0 * u
    @test @inferred(etp(2018)) === 2.0f0 * u
    @test @inferred(etp(2019)) === 1.0f0 * u

    # Unitful vector knots
    t = [0.0, 1.0] * u"m"
    itp = linear_interpolation(t, y)
    @test itp(0.5u"m") == 1.5f0 * u
    @test itp([0.5, 0.75]u"m") == [1.5f0, 1.75f0] * u

    # issues #370, #424
    itp = linear_interpolation([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], extrapolation_bc=missing)
    @test ismissing(itp(4.0))
    @test @inferred(eltype(itp),itp(3.0)) === 3.0
    @test @inferred(eltype(itp),itp(4.0)) === missing
end
