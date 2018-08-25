@testset "NoInterp" begin
    a = reshape(1:12, 3, 4)
    ai = interpolate(a, NoInterp())
    @test eltype(ai) == Int
    check_axes(ai, a)
    check_inbounds_values(ai, a)
    check_oob(ai)
    @test_throws InexactError ai(2.2, 2)
    @test_throws InexactError ai(2, 2.2)

    # ae = extrapolate(ai, NaN)
    # @test eltype(ae) == Float64
    # @test ae[1,1] === 1.0
    # @test ae[0,1] === NaN
    # @test_throws InexactError ae(1.5,2)
end
