using Interpolations, DualNumbers, Test, LinearAlgebra

@testset "ScalingDimspecTests" begin
    xs = -pi:(2pi/10):pi-2pi/10
    ys = -2:.1:2
    f(x,y) = sin(x) * y^2

    itp = interpolate(Float64[f(x,y) for x in xs, y in ys], (BSpline(Quadratic(Periodic(OnGrid()))), BSpline(Linear())))
    sitp = scale(itp, xs, ys)

    for (ix,x) in enumerate(xs), (iy,y) in enumerate(ys)
        @test ≈(sitp(x,y),f(x,y),atol=sqrt(eps(1.0)))

        g = Interpolations.gradient(sitp, x, y)
        fx = epsilon(sitp(dual(x,1), dual(y,0)))
        fy = epsilon(sitp(dual(x,0), dual(y,1)))

        @test ≈(g[1],fx,atol=sqrt(eps(1.0)))
        @test ≈(g[2],fy,atol=sqrt(eps(1.0)))
    end
end
