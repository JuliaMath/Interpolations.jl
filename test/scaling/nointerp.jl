module ScalingNoInterpTests

using Interpolations, Base.Test

xs = -pi:2pi/10:pi
f1(x) = sin(x)
f2(x) = cos(x)
f3(x) = sin(x) .* cos(x)
f(x,y) = y == 1 ? f1(x) : (y == 2 ? f2(x) : (y == 3 ? f3(x) : error("invalid value for y (must be 1, 2 or 3, you used $y)")))
ys = 1:3

A = hcat(map(f1, xs), map(f2, xs), map(f3, xs))

itp = interpolate(A, (BSpline(Quadratic(Periodic())), NoInterp()), OnGrid())
sitp = scale(itp, xs, ys)

for (ix,x0) in enumerate(xs[1:end-1]), y0 in ys
    x,y = x0, y0
    @test ≈(sitp[x,y],f(x,y),atol=0.05)
end

@test length(gradient(sitp, pi/3, 2)) == 1


end
