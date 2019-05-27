module TypeInstantiationTests

using Interpolations, Test

# NO DIMSPECS
# tests that we forward types correctly to the instance constructors
const A1 = rand(15)
const A2 = rand(15, 10)

# B-splines
for GT in (OnCell, OnGrid)
    # ...without boundary conditions
    for D in (Constant, Linear)
        @inferred(interpolate(A1, BSpline(D), GT))
        @inferred(interpolate(A1, BSpline(D), GT()))
        @inferred(interpolate(A1, BSpline(D()), GT))
        @inferred(interpolate(A2, BSpline(D), GT))
        @inferred(interpolate(A2, BSpline(D), GT()))
        @inferred(interpolate(A2, BSpline(D()), GT))

        @inferred(interpolate!(copy(A1), BSpline(D), GT))
        @inferred(interpolate!(copy(A1), BSpline(D), GT()))
        @inferred(interpolate!(copy(A1), BSpline(D()), GT))
        @inferred(interpolate!(copy(A2), BSpline(D), GT))
        @inferred(interpolate!(copy(A2), BSpline(D), GT()))
        @inferred(interpolate!(copy(A2), BSpline(D()), GT))
    end

    # Quadratic
    for BC in (Flat, Line, Periodic, Reflect, Free)
        if BC != InPlace
            @inferred(interpolate(A1, BSpline(Quadratic(BC)), GT))
            @inferred(interpolate(A1, BSpline(Quadratic(BC())), GT))
            @inferred(interpolate(A1, BSpline(Quadratic(BC)), GT()))
            @inferred(interpolate(A2, BSpline(Quadratic(BC)), GT))
            @inferred(interpolate(A2, BSpline(Quadratic(BC())), GT))
            @inferred(interpolate(A2, BSpline(Quadratic(BC)), GT()))
        end

        @inferred(interpolate!(copy(A1), BSpline(Quadratic(BC)), GT))
        @inferred(interpolate!(copy(A1), BSpline(Quadratic(BC())), GT))
        @inferred(interpolate!(copy(A1), BSpline(Quadratic(BC)), GT()))
        @inferred(interpolate!(copy(A2), BSpline(Quadratic(BC)), GT))
        @inferred(interpolate!(copy(A2), BSpline(Quadratic(BC())), GT))
        @inferred(interpolate!(copy(A2), BSpline(Quadratic(BC)), GT()))
    end
end

# Gridded
const knots1 = (sort(rand(15)),)
const knots2 = (sort(rand(15)), sort(rand(10)))
for D in (Constant, Linear)
    @inferred(interpolate(knots1, A1, Gridded(D)))
    @inferred(interpolate(knots2, A2, Gridded(D)))

    @inferred(interpolate!(knots1, copy(A1), Gridded(D)))
    @inferred(interpolate!(knots2, copy(A2), Gridded(D)))
end

# DIMSPECS

# test that constructing dimspecs work
for T in (
    Tuple{OnGrid,OnCell},
    Tuple{BSpline{Linear},BSpline{Quadratic{Free}}}
    )
    @test isa(@inferred(Interpolations.construct_instance(T)), T)
end
# sample one dimspec for each interpolation constructor to see that it
# calls the construct_instance() function correctly
@inferred(interpolate(A2, Tuple{BSpline{Linear},BSpline{Constant}}, Tuple{OnGrid,OnCell}))
@inferred(interpolate!(copy(A2), Tuple{BSpline{Linear},BSpline{Constant}}, Tuple{OnGrid,OnCell}))
@inferred(interpolate(A2, Tuple{BSpline{Linear},BSpline{Constant}}, (OnGrid(),OnCell())))
@inferred(interpolate!(copy(A2), Tuple{BSpline{Linear},BSpline{Constant}}, (OnGrid(),OnCell())))
@inferred(interpolate(A2, (BSpline(Linear()),BSpline(Constant())), Tuple{OnGrid,OnCell}))
@inferred(interpolate!(copy(A2), (BSpline(Linear()),BSpline(Constant())), Tuple{OnGrid,OnCell}))

end
