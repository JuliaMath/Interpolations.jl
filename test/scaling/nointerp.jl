using Interpolations, Test, LinearAlgebra, Random

@testset "ScalingNoInterpTests" begin
    xs = -pi:2pi/10:pi
    f1(x) = sin(x)
    f2(x) = cos(x)
    f3(x) = sin(x) * cos(x)
    f(x,y) = y == 1 ? f1(x) : (y == 2 ? f2(x) : (y == 3 ? f3(x) : error("invalid value for y (must be 1, 2 or 3, you used $y)")))
    ys = 1:3

    A = hcat(map(f1, xs), map(f2, xs), map(f3, xs))

    itp = interpolate(A, (BSpline(Quadratic(Periodic(OnGrid()))), NoInterp()))
    sitp = scale(itp, xs, ys)

    for (ix,x0) in enumerate(xs[1:end-1]), y0 in ys
        x,y = x0, y0
        @test sitp(x,y) ≈ f(x,y) atol=0.05
    end

    @test length(Interpolations.gradient(sitp, pi/3, 2)) == 1

    xs = range(-pi, stop=pi, length=60)[1:end-1]
    A = hcat(map(f1, xs), map(f2, xs), map(f3, xs))
    itp = interpolate(A, (BSpline(Cubic(Periodic(OnGrid()))), NoInterp()))
    sitp = scale(itp, xs, ys)

    for x in xs, y in ys
        if y in (1,2)
            h = @inferred(Interpolations.hessian(sitp, x, y))
            @test h[1] ≈ -f(x, y) atol=0.01
        else # y==3
            h = @inferred(Interpolations.hessian(sitp, x, y))
            @test h[1] ≈ -4*f(x, y) atol=0.01
        end
    end

    @test length(Interpolations.hessian(sitp, pi/3, 2)) == 1

    A = hcat(map(f1, xs), map(f2, xs), map(f3, xs))
    itp = interpolate(A', (NoInterp(), BSpline(Cubic(Periodic(OnGrid())))))
    sitp = scale(itp, ys, xs)
    h(y, x) = Interpolations.hessian(sitp, y, x)
    h(1, xs[10])
    for x in xs[6:end-6], y in ys
        if y in (1,2)
            @test h(y, x)[1] ≈ -f(x, y) atol=0.05
        elseif y==3
            @test h(y, x)[1] ≈ -4*f(x, y) atol=0.05
        end
    end

    # check for case where initial/middle indices are NoInterp but later ones are <:BSpline
    isdefined(Random, :seed!) ? Random.seed!(1234) : srand(1234) # `srand` was renamed to `seed!`
    z0 = rand(10,10)
    za = copy(z0)
    zb = copy(z0')

    itpa = interpolate(za, (BSpline(Linear()), NoInterp()))
    itpb = interpolate(zb, (NoInterp(), BSpline(Linear())))

    rng = range(1.0, stop=19.0, length=10)
    sitpa = scale(itpa, rng, 1:10)
    sitpb = scale(itpb, 1:10, rng)
    @test Interpolations.gradient(sitpa, 3.0, 3) ==  Interpolations.gradient(sitpb, 3, 3.0)


end
