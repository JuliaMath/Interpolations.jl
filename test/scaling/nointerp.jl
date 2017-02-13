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

if VERSION < v"0.5.0-dev" # Test.with_handler was removed in 0.5

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
    Test.with_handler(message_is("Must index into 2-dimensional scaled interpolation object with exactly 2 indices (you used 1)")) do
        @test gradient(sitp, 2.3)
    end
    Test.with_handler(message_is("The length of the provided gradient vector (2) did not match the number of interpolating dimensions (1)")) do
        @test gradient!(Array{Float64}( 2), sitp, 2.3, 2)
    end

end

end
