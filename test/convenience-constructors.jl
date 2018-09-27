module ConvenienceConstructorTests

using Interpolations
using Test
using Base.Cartesian

# unit test setup
XMIN = 2
XMAX = 10
YMIN = 1
YMAX = 8
ΔX = .1
ΔY = .5
XLEN = convert(Integer, floor((XMAX - XMIN)/ΔX) + 1)
YLEN = convert(Integer, floor((YMAX - YMIN)/ΔY) + 1)

@testset "1d-interpolations" begin
    @testset "1d-regular-grids" begin
        xs = XMIN:ΔX:XMAX
        f(x) = log(x)
        A = [f(x) for x in xs]
        interp = LinearInterpolation(xs, A) # using convenience constructor
        interp_full = extrapolate(scale(interpolate(A, BSpline(Linear())), xs), Throw()) # using full constructor

        @test typeof(interp) == typeof(interp_full)
        @test interp(XMIN) ≈ f(XMIN)
        @test interp(XMAX) ≈ f(XMAX)
        @test interp(XMIN + ΔX) ≈ f(XMIN + ΔX)
        @test interp(XMAX - ΔX) ≈ f(XMAX - ΔX)
        @test interp(XMIN + ΔX / 2) ≈ f(XMIN + ΔX / 2) atol=.1
        @test_throws BoundsError interp(XMIN - ΔX / 2)
        @test_throws BoundsError interp(XMAX + ΔX / 2)
    end

    @testset "1d-regular-grids-cubic" begin
        xs = XMIN:ΔX:XMAX
        f(x) = log(x)
        A = [f(x) for x in xs]
        interp = CubicSplineInterpolation(xs, A)
        interp_full = extrapolate(scale(interpolate(A, BSpline(Cubic(Line(OnGrid())))), xs), Throw())

        @test typeof(interp) == typeof(interp_full)
        @test interp(XMIN) ≈ f(XMIN)
        @test interp(XMAX) ≈ f(XMAX)
        @test interp(XMIN + ΔX) ≈ f(XMIN + ΔX)
        @test interp(XMAX - ΔX) ≈ f(XMAX - ΔX)
        @test interp(XMIN + ΔX / 2) ≈ f(XMIN + ΔX / 2) atol=.1
        @test_throws BoundsError interp(XMIN - ΔX / 2)
        @test_throws BoundsError interp(XMAX + ΔX / 2)
    end

    @testset "1d-irregular-grids" begin
        xs = [x^2 for x in XMIN:ΔX:XMAX]
        xmin = xs[1]
        xmax = xs[XLEN]
        f(x) = log(x)
        A = [f(x) for x in xs]
        interp = LinearInterpolation(xs, A)
        interp_full = extrapolate(interpolate((xs, ), A, Gridded(Linear())), Throw())

        @test typeof(interp) == typeof(interp_full)
        @test interp(xmin) ≈ f(xmin)
        @test interp(xmax) ≈ f(xmax)
        @test interp(xs[2]) ≈ f(xs[2])
        @test interp(xmin + ΔX / 2) ≈ f(xmin + ΔX / 2) atol=.1
        @test_throws BoundsError interp(xmin - ΔX / 2)
        @test_throws BoundsError interp(xmax + ΔX / 2)
    end

    @testset "1d-handling-extrapolation" begin
        xs = XMIN:ΔX:XMAX
        f(x) = log(x)
        A = [f(x) for x in xs]
        ΔA_l = A[2] - A[1]
        ΔA_h = A[end] - A[end - 1]
        x_lower = XMIN - ΔX
        x_higher = XMAX + ΔX

        extrap = LinearInterpolation(xs, A, extrapolation_bc = Line())
        extrap_full = extrapolate(scale(interpolate(A, BSpline(Linear())), xs), Line())

        @test typeof(extrap) == typeof(extrap_full)
        @test extrap(x_lower) ≈ A[1] - ΔA_l
        @test extrap(x_higher) ≈ A[end] + ΔA_h
    end
end

@testset "2d-interpolations" begin
    @testset "2d-regular-grids" begin
        xs = XMIN:ΔX:XMAX
        ys = YMIN:ΔY:YMAX
        f(x, y) = log(x+y)
        A = [f(x,y) for x in xs, y in ys]
        interp = LinearInterpolation((xs, ys), A)
        interp_full = extrapolate(scale(interpolate(A, BSpline(Linear())), xs, ys), Throw())

        @test typeof(interp) == typeof(interp_full)
        @test interp(XMIN,YMIN) ≈ f(XMIN,YMIN)
        @test interp(XMIN,YMAX) ≈ f(XMIN,YMAX)
        @test interp(XMAX,YMIN) ≈ f(XMAX,YMIN)
        @test interp(XMAX,YMAX) ≈ f(XMAX,YMAX)
        @test interp(XMIN + ΔX,YMIN) ≈ f(XMIN + ΔX,YMIN)
        @test interp(XMIN,YMIN + ΔY) ≈ f(XMIN,YMIN + ΔY)
        @test interp(XMIN + ΔX,YMIN + ΔY) ≈ f(XMIN + ΔX,YMIN + ΔY)
        @test interp(XMIN + ΔX / 2,YMIN + ΔY / 2) ≈ f(XMIN + ΔX / 2,YMIN + ΔY / 2) atol=.1
        @test_throws BoundsError interp(XMIN - ΔX / 2,YMIN - ΔY / 2)
        @test_throws BoundsError interp(XMIN - ΔX / 2,YMIN + ΔY / 2)
        @test_throws BoundsError interp(XMIN + ΔX / 2,YMIN - ΔY / 2)
        @test_throws BoundsError interp(XMAX + ΔX / 2,YMAX + ΔY / 2)
    end

    @testset "2d-regular-grids-cubic" begin
        xs = XMIN:ΔX:XMAX
        ys = YMIN:ΔY:YMAX
        f(x, y) = log(x+y)
        A = [f(x,y) for x in xs, y in ys]
        interp = CubicSplineInterpolation((xs, ys), A)
        interp_full = extrapolate(scale(interpolate(A, BSpline(Cubic(Line(OnGrid())))), xs, ys), Throw())

        @test typeof(interp) == typeof(interp_full)
        @test interp(XMIN,YMIN) ≈ f(XMIN,YMIN)
        @test interp(XMIN,YMAX) ≈ f(XMIN,YMAX)
        @test interp(XMAX,YMIN) ≈ f(XMAX,YMIN)
        @test interp(XMAX,YMAX) ≈ f(XMAX,YMAX)
        @test interp(XMIN + ΔX,YMIN) ≈ f(XMIN + ΔX,YMIN)
        @test interp(XMIN,YMIN + ΔY) ≈ f(XMIN,YMIN + ΔY)
        @test interp(XMIN + ΔX,YMIN + ΔY) ≈ f(XMIN + ΔX,YMIN + ΔY)
        @test interp(XMIN + ΔX / 2,YMIN + ΔY / 2) ≈ f(XMIN + ΔX / 2,YMIN + ΔY / 2) atol=.1
        @test_throws BoundsError interp(XMIN - ΔX / 2,YMIN - ΔY / 2)
        @test_throws BoundsError interp(XMIN - ΔX / 2,YMIN + ΔY / 2)
        @test_throws BoundsError interp(XMIN + ΔX / 2,YMIN - ΔY / 2)
        @test_throws BoundsError interp(XMAX + ΔX / 2,YMAX + ΔY / 2)
    end

    @testset "2d-irregular-grids" begin
        xs = [x^2 for x in XMIN:ΔX:XMAX]
        ys = [y^2 for y in YMIN:ΔY:YMAX]
        xmin = xs[1]
        xmax = xs[XLEN]
        ymin = ys[1]
        ymax = ys[YLEN]
        f(x, y) = log(x+y)
        A = [f(x,y) for x in xs, y in ys]
        interp = LinearInterpolation((xs, ys), A)
        interp_full = extrapolate(interpolate((xs, ys), A, Gridded(Linear())), Throw())

        @test typeof(interp) == typeof(interp_full)
        @test interp(xmin,ymin) ≈ f(xmin,ymin)
        @test interp(xmin,ymax) ≈ f(xmin,ymax)
        @test interp(xmax,ymin) ≈ f(xmax,ymin)
        @test interp(xmax,ymax) ≈ f(xmax,ymax)
        @test interp(xs[2],ymin) ≈ f(xs[2],ymin)
        @test interp(xmin,ys[2]) ≈ f(xmin,ys[2])
        @test interp(xs[2],ys[2]) ≈ f(xs[2],ys[2])
        @test interp(xmin + ΔX / 2,ymin + ΔY / 2) ≈ f(xmin + ΔX / 2,ymin + ΔY / 2) atol=.1
        @test_throws BoundsError interp(xmin - ΔX / 2,ymin - ΔY / 2)
        @test_throws BoundsError interp(xmin - ΔX / 2,ymin + ΔY / 2)
        @test_throws BoundsError interp(xmin + ΔX / 2,ymin - ΔY / 2)
        @test_throws BoundsError interp(xmax + ΔX / 2,ymax + ΔY / 2)
    end

    @testset "2d-handling-extrapolation" begin
        xs = XMIN:ΔX:XMAX
        ys = YMIN:ΔY:YMAX
        f(x, y) = log(x+y)
        A = [f(x,y) for x in xs, y in ys]
        ΔA_l = A[2, 1] - A[1, 1]
        ΔA_h = A[end, end] - A[end - 1, end]
        x_lower = XMIN - ΔX
        x_higher = XMAX + ΔX
        y_lower = YMIN - ΔY
        y_higher = YMAX + ΔY

        extrap = LinearInterpolation((xs, ys), A, extrapolation_bc = (Line(), Flat()))
        extrap_full = extrapolate(scale(interpolate(A, BSpline(Linear())), xs, ys), (Line(), Flat()))

        @test typeof(extrap) == typeof(extrap_full)
        @test extrap(x_lower, y_lower) ≈ A[1, 1] - ΔA_l
        @test extrap(x_higher, y_higher) ≈ A[end, end] + ΔA_h
    end

    @testset "issue #230" begin # at least, I think this is what issue #230 is really about
        f(x,y) = log(x+y)
        xs = 1:5
        ys = 2:0.1:5
        A = [f(x,y) for x in xs, y in ys]
        itp = LinearInterpolation((xs, ys), A)
        for (i, j) in zip(Iterators.product(xs, ys), eachindex(A))
            @test itp(i...) ≈ A[j]
        end
    end
end

end
