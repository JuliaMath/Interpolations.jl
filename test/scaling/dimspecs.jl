module ScalingDimspecTests

using Interpolations, DualNumbers, Base.Test

xs = -pi:(2pi/10):pi-2pi/10
ys = -2:.1:2
f(x,y) = sin(x) * y^2

itp = interpolate(Float64[f(x,y) for x in xs, y in ys], (BSpline(Quadratic(Periodic())), BSpline(Linear())), OnGrid())
sitp = scale(itp, xs, ys)

for (ix,x) in enumerate(xs), (iy,y) in enumerate(ys)
    @test_approx_eq_eps sitp[x,y] f(x,y) sqrt(eps(1.0))

    g = gradient(sitp, x, y)
    fx = epsilon(sitp[dual(x,1), dual(y,0)])
    fy = epsilon(sitp[dual(x,0), dual(y,1)])

    @test_approx_eq_eps g[1] fx sqrt(eps(1.0))
    @test_approx_eq_eps g[2] fy sqrt(eps(1.0))
end

end