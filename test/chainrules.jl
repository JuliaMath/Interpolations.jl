using Zygote: Zygote
@testset "ChainRulesCore" begin
    # 1D example
    x = 1:10
    y = sin.(x)

    # Simple test that rrules give equivalent results to internal gradient function

    itp = interpolate(y,BSpline(Linear()))

    @test Zygote.gradient(itp, 1) == Tuple(Interpolations.gradient(itp, 1))

    # 2D example
    x2 = [i*j for i=1:5, j=1:5]
    y2 = sin.(x2)

    itp2 = interpolate(y2, BSpline(Cubic(Line(OnGrid()))))

    @test Zygote.gradient(itp2,1,2) == Tuple(Interpolations.gradient(itp2,1,2))

    # Matrix-valued
    y=[[[0 0; 0 0]] [[1 2; 0 0]];
       [[0 0; 3 4]] [[1 2; 3 4]]]
    itp_m = interpolate(y, BSpline(Linear()))
    @test Zygote.jacobian(itp_m, 1,1) == Tuple((x->[x...]).(Interpolations.gradient(itp_m, 1,1)))

end
