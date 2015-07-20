module ScalingDimspecTests

using Interpolations, DualNumbers, Base.Test

xs = -pi:(2pi/10):pi-2pi/10
ys = -2:.1:2
f(x,y) = sin(x) * y^2

itp = interpolate(Float64[f(x,y) for x in xs, y in ys], Tuple{BSpline(Quadratic(Periodic)), BSpline(Linear)}, OnGrid)
sitp = scale(itp, xs, ys)

# Don't test too near the edges until issue #64 is resolved
for (ix,x0) in enumerate(xs[5:end-5]), (iy,y0) in enumerate(ys[2:end-1])
    x, y = x0 + 2pi/20, y0 + .05
    @test_approx_eq sitp[x0, y0] f(x0,y0)
    @test_approx_eq_eps sitp[x0, y0] f(x0,y0) 0.05

    g = gradient(sitp, x, y)
    fx = epsilon(f(dual(x,1), y))
    fy = (f(x, ys[iy+2]) - f(x, ys[iy+1])) / (ys[iy+2] - ys[iy+1])

    @test_approx_eq_eps g[1] fx 0.15
    @test_approx_eq_eps g[2] fy 0.05 # gradients for linear interpolation is "easy"
end

end
