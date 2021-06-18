
@testset "Periodic extrapolation" begin
    xmax = 7
    f(x) = sin((x-3)*2pi/xmax - 1)
    A = Float64[f(x) for x in 1:xmax] # Does not include the periodic image

    itp0 = interpolate(A, BSpline(Constant(Periodic())))
    itp1 = interpolate(A, BSpline(Linear(Periodic())))
    itp2 = interpolate(A, BSpline(Quadratic(Periodic(OnCell()))))
    itp3 = interpolate(A, BSpline(Cubic(Periodic(OnCell()))))

    etp0 = extrapolate(itp0, Periodic())
    etp1 = extrapolate(itp1, Periodic())
    etp2 = extrapolate(itp2, Periodic())
    etp3 = extrapolate(itp3, Periodic())

    for x in -xmax:.4:2*xmax
        @test etp0(x) ≈ f(x) atol=0.5
        @test etp1(x) ≈ f(x) atol=0.1
        @test etp2(x) ≈ f(x) atol=0.01
        @test etp3(x) ≈ f(x) atol=0.003
    end
end