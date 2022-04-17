using Zygote
@testset "ChainRulesCore" begin
    # 1D example
    x = 1:10
    y = sin.(x)

    # Simple test that rrules give equivalent results to internal gradient function

    itp = interpolate(y,BSpline(Linear()))

    if VERSION ≥ v"1.3"
        @test Zygote.gradient(itp, 1) == Tuple(Interpolations.gradient(itp, 1))
    else
        @test_skip Zygote.gradient(itp, 1) == Tuple(Interpolations.gradient(itp, 1))
    end

    # 2D example
    x2 = [i*j for i=1:5, j=1:5]
    y2 = sin.(x2)

    itp2 = interpolate(y2, BSpline(Cubic(Line(OnGrid()))))

    if VERSION ≥ v"1.3"
        @test Zygote.gradient(itp2,1,2) == Tuple(Interpolations.gradient(itp2,1,2))
    else
        @test_skip Zygote.gradient(itp2,1,2) == Tuple(Interpolations.gradient(itp2,1,2))
    end

    # Matrix-valued
    y=[[[0 0; 0 0]]; [[1 2; 0 0]];;
       [[0 0; 3 4]]; [[1 2; 3 4]]]
    itp_m = interpolate(y, BSpline(Linear()))
    if VERSION ≥ v"1.3"
        @test Zygote.jacobian(itp_m, 1,1) == Tuple((x->[x...]).(Interpolations.gradient(itp_m, 1,1)))
    else
        @test_skip Zygote.jacobian(itp_m, 1,1) == Tuple((x->[x...]).(Interpolations.gradient(itp_m, 1,1)))
    end

end
