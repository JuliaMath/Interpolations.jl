using Test
using Interpolations

@testset "Pipes" begin
    for BC in (Throw, Flat, Line, Free, Reflect, InPlace, InPlaceQ, Periodic)
        for GT in (OnGrid, OnCell)
            @test BC(GT) == BC(GT())
            @test GT |> BC == BC(GT())
            for D in (Cubic, Quadratic)
                @test GT |> BC |> D == D(BC(GT()))
                @test GT |> BC |> D{BC{GT}} == D(BC(GT()))
            end
        end
    end
    @test BSpline() == BSpline(Linear())
end
