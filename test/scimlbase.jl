using Interpolations, SciMLBase, StaticArrays

@testset "SciMLBase should_warn_paramtype" begin
    Xs = 10.0:2.0:20.0
    # want Y-values that are
    #   - Concrete
    Ys_float = Xs .^ 2
    #   - Abstract
    Ys_number = Vector{Real}(Ys_float)
    #   - Array of concrete
    Ys_arr = map(x -> SVector(x^2, 2 * x - 1), Xs)
    #   - Array of abstract
    Ys_arr_number = map(((i, x),) -> SVector{2, Real}(i, x), enumerate(Xs))
    bspline_concrete_scalar = interpolate(Ys_float, BSpline(Linear()))
    bspline_abstract_scalar = interpolate(Ys_number, BSpline(Linear()))
    bspline_concrete_array = interpolate(Ys_arr, BSpline(Linear()))
    bspline_abstract_array = interpolate(Ys_arr_number, BSpline(Linear()))

    gridded_concrete_scalar = interpolate((Xs,), Ys_float, Gridded(Linear()))
    gridded_abstract_scalar = interpolate((Xs,), Ys_number, Gridded(Linear()))
    gridded_concrete_array = interpolate((Xs,), Ys_arr, Gridded(Linear()))
    gridded_abstract_array = interpolate((Xs,), Ys_arr_number, Gridded(Linear()))

    # should not warn on concrete number
    @test SciMLBase.should_warn_paramtype(bspline_concrete_scalar) == false
    @test SciMLBase.should_warn_paramtype(gridded_concrete_scalar) == false
    # should not warn on concrete array eltype
    @test SciMLBase.should_warn_paramtype(bspline_concrete_array) == false
    @test SciMLBase.should_warn_paramtype(gridded_concrete_array) == false
    # should warn on abstract number
    @test SciMLBase.should_warn_paramtype(bspline_abstract_scalar)
    @test SciMLBase.should_warn_paramtype(gridded_abstract_scalar)
    # should warn on abstract array eltype
    @test SciMLBase.should_warn_paramtype(bspline_abstract_array)
    @test SciMLBase.should_warn_paramtype(gridded_abstract_array)
end
