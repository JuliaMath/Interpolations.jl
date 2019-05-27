using Test, Interpolations, OffsetArrays

@testset "ExtrapNon1" begin

    f(x) = sin((x-3)*2pi/9 - 1)
    xinds = -3:6
    A = OffsetArray(Float64[f(x) for x in xinds], xinds)
    itpg = interpolate(A, BSpline(Linear()))

    schemes = (
        Flat,
        Line,
        Reflect,
        Periodic
    )

    for etp in map(E -> extrapolate(itpg, E()), schemes), x in xinds
        @test parent(etp) === itpg
        @test @inferred(getindex(etp, x)) ≈ A[x]
    end

    g(y) = (y/100)^3
    yinds = 2:5
    A = OffsetArray(Float64[f(x)*g(y) for x in xinds, y in yinds], xinds, yinds)
    itp2 = interpolate(A, BSpline(Linear()))

    for (etp2,E) in map(E -> (extrapolate(itp2, E()), E), schemes)
        @test parent(etp2) === itp2
        E == Periodic && continue  # g isn't periodic
        for y in yinds, x in xinds
            @test @inferred(getindex(etp2, x, y)) ≈ A[x, y]
        end
    end

end
