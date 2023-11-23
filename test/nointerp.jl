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

@testset "Stability of mixtrue with NoInterp and Interp" begin
    A = zeros(Float64, 5, 5, 5, 5, 5, 5, 5)
    st = BSpline(Quadratic(Reflect(OnCell()))), NoInterp(),
         BSpline(Linear()), NoInterp(),
         BSpline(Quadratic()), NoInterp(),
         BSpline(Quadratic(Reflect(OnCell())))
    itp = interpolate(A, st)
    @test (@inferred Interpolations.hessian(itp, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)) == zeros(4,4)
    @test (@inferred Interpolations.gradient(itp, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)) == zeros(4)
    @test (@inferred itp(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)) == 0
end
