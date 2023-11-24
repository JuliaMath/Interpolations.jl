using Test
using Interpolations

@testset "Pipes" begin
    for BC in (Throw, Flat, Line, Free, Reflect, InPlace, InPlaceQ, Periodic)
        for GT in (OnGrid, OnCell)
            @test BC(GT) == BC(GT())
            @test BC{GT}(GT) == BC{GT}(GT())
            @test GT |> BC == BC(GT())
            for D in (Cubic, Quadratic)
                @test GT |> BC |> D == D(BC(GT()))
                @test GT |> BC |> D{BC{GT}} == D(BC(GT()))
                @test D{BC}(BC) == D{BC}(BC())
            end
        end
        @test Linear(Throw{OnGrid}) == Linear(Throw{OnGrid}())
    end
    @test BSpline() == BSpline(Linear())
    @test_throws ArgumentError Cubic{Throw}(Flat)
    @test_throws ArgumentError Quadratic{Throw}(Flat)
end
