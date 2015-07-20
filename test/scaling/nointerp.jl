module ScalingNoInterpTests

using Interpolations, Base.Test

xs = -pi:2pi/10:pi
f1(x) = sin(x)
f2(x) = cos(x)
f3(x) = sin(x) .* cos(x)
f(x,y) = y == 1 ? f1(x) : (y == 2 ? f2(x) : (y == 3 ? f3(x) : error("invalid value for y (must be 1, 2 or 3, you used $y)")))
ys = 1:3

A = hcat(f1(xs), f2(xs), f3(xs))

itp = interpolate(A, Tuple{BSpline(Quadratic(Periodic)), NoInterp}, OnGrid)
sitp = scale(itp, xs, ys)

for (ix,x0) in enumerate(xs[1:end-1]), y0 in ys
    x,y = x0, y0
    @test_approx_eq_eps sitp[x,y] f(x,y) .05
end

# Test error messages for incorrect initialization
function message_is(message)
	r -> r.err.msg == message || error("Incorrect error message: expected '$message' but was '$(r.err.msg)'")
end
Test.with_handler(message_is("Must scale 2-dimensional interpolation object with exactly 2 ranges (you used 1)")) do
	@test scale(itp, xs)
end
Test.with_handler(message_is("NoInterp dimension 2 must be scaled with unit range 1:3")) do
	@test scale(itp, xs, -1:1)
end
Test.with_handler(message_is("The length of the range in dimension 1 (8) did not equal the size of the interpolation object in that direction (11)")) do
	@test scale(itp, -pi:2pi/7:pi, 1:3)
end
Test.with_handler(message_is("Must index into 2-dimensional scaled interpolation object with exactly 2 indices (you used 1)")) do
	@test sitp[2.3]
end
end
