using Random
using Test

function cubic_hermite_test()
    # Need at least two samples:
    xs = ones(1)
    ys = ones(1)
    dydxs = ones(1)
    @test_throws DomainError CubicHermite(xs, ys, dydxs)

    # All vectors must be the same size:
    xs = ones(2)
    xs[2] = 2
    ys = ones(2)
    dydxs = ones(3)
    @test_throws DomainError CubicHermite(xs, ys, dydxs)

    xs = ones(5)
    for i = 2:length(xs)
        xs[i] = i
    end
    ys = ones(5)
    dydxs = zeros(5)
    ch = CubicHermite(xs, ys, dydxs)
    for i = 1:length(xs)-1
        x = Float64(i)
        @test ch(x) == 1.0
        @test Interpolations.gradient(ch, x) == 0.0
        @test Interpolations.hessian(ch, x) == 0.0
        x = x + 0.25
        @test ch(x) == 1.0
        @test Interpolations.gradient(ch, x) == 0.0
        @test Interpolations.hessian(ch, x) == 0.0
        x = x + 0.25
        @test ch(x) == 1.0
        @test Interpolations.gradient(ch, x) == 0.0
        @test Interpolations.hessian(ch, x) == 0.0
        x = x + 0.25
        @test ch(x) == 1.0
        @test Interpolations.gradient(ch, x) == 0.0
        @test Interpolations.hessian(ch, x) == 0.0
    end
    # Ensure that the ight endpoint query doesn't read past the end of the array:
    @test ch(5.0) == 1.0
    @test Interpolations.gradient(ch, 5.0) == 0.0
    @test Interpolations.hessian(ch, 5.0) == 0.0

    # Now linear functions:
    a = 7.2
    b = 9.6
    for i = 1:length(xs)
        ys[i] = a * xs[i] + b
        dydxs[i] = a
    end
    ch = CubicHermite(xs, ys, dydxs)
    for i = 1:length(xs)-1
        x = Float64(i)
        @test ch(x) ≈ a * x + b
        @test Interpolations.gradient(ch, x) ≈ a
        @test abs(Interpolations.hessian(ch, x)) < 3e-14
        x = x + 0.25
        @test ch(x) ≈ a * x + b
        @test Interpolations.gradient(ch, x) ≈ a
        @test abs(Interpolations.hessian(ch, x)) < 3e-14
        x = x + 0.25
        @test ch(x) ≈ a * x + b
        @test Interpolations.gradient(ch, x) ≈ a
        @test abs(Interpolations.hessian(ch, x)) < 3e-14
        x = x + 0.25
        @test ch(x) ≈ a * x + b
        @test Interpolations.gradient(ch, x) ≈ a
        @test abs(Interpolations.hessian(ch, x)) < 3e-14
    end
    @test ch(last(xs)) ≈ a * last(xs) + b
    @test Interpolations.gradient(ch, last(xs)) ≈ a
    @test abs(Interpolations.hessian(ch, last(xs))) < 3e-14

    # Now the interpolation condition:
    xs = zeros(50)
    ys = zeros(50)
    dydxs = zeros(50)
    xs[1] = 0.0
    ys[1] = randn()
    dydxs[1] = randn()
    for i = 2:50
        xs[i] = xs[i-1] + rand() + 0.1
        ys[i] = randn()
        dydxs[i] = randn()
    end

    ch = CubicHermite(xs, ys, dydxs)

    for i = 1:50
        @test ch(xs[i]) ≈ ys[i]
        @test Interpolations.gradient(ch, xs[i]) ≈ dydxs[i]
    end

    # Now quadratics:
    a = 1.0 / 8
    b = -3
    c = -2
    for i = 1:50
        ys[i] = a * xs[i] * xs[i] + b * xs[i] + c
        dydxs[i] = 2 * a * xs[i] + b
    end
    ch = CubicHermite(xs, ys, dydxs)
    for i = 1:200
        x = rand() * last(xs)
        @test ch(x) ≈ a * x * x + b * x + c
        @test Interpolations.gradient(ch, x) ≈ 2 * a * x + b
        @test Interpolations.hessian(ch, x) ≈ 2 * a
    end
    x = last(xs)
    @test ch(x) ≈ a * x * x + b * x + c
    @test Interpolations.gradient(ch, x) ≈ 2 * x * a + b
    @test Interpolations.hessian(ch, x) ≈ 2 * a
    # Cannot extrapolate:
    @test_throws DomainError ch(x + 0.1)
    @test_throws DomainError ch(xs[1] - 0.1)
end

@testset "Hermite" begin
    cubic_hermite_test()
end
