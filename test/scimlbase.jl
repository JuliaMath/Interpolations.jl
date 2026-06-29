using SciMLBase, StaticArrays, Interpolations

@testset "SciMLBase should_warn_paramtype" begin
    # Abstract range
    Xs_float = 0.0:2.0:10.0
    Ys_float = Xs_float .^ 2
    # this might be clunky
    Ys_number = Vector{Number}(Ys_float)
    Ys_arr = map(x -> SVector(x^2, 2 * x - 1), Xs_float)
    Ys_arr_number = map(((i, x),) -> SVector{2, Number}(i, x), enumerate(Xs_float))
    i_float = interpolate(Ys_float, BSpline(Linear()))
    i_number = interpolate(Ys_number, BSpline(Linear()))
    i_arr_float = interpolate(Ys_arr, BSpline(Linear()))
    i_arr_number = interpolate(Ys_arr_number, NoInterp())
    # should not warn on concrete number
    @test SciMLBase.should_warn_paramtype(i_float) == false
    # should not warn on concrete array eltype
    @test SciMLBase.should_warn_paramtype(i_arr_float) == false
    # call operator on abstract number interpolation seems to return f64
    @test SciMLBase.should_warn_paramtype(i_number) == false
    # should warn on abstract array eltype
    @test SciMLBase.should_warn_paramtype(i_arr_number)
end
