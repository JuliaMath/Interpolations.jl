using DualNumbers
using Interpolations
using Test

@testset "Extrapolation" begin

    f(x) = sin((x-3)*2pi/9 - 1)
    xmax = 10
    A = Float64[f(x) for x in 1:xmax]

    itpg = interpolate(A, BSpline(Linear()))

    etpg = extrapolate(itpg, Flat())
    @test typeof(etpg) <: AbstractExtrapolation

    @test etpg(-3) == etpg(-4.5) == etpg(0.9) == etpg(1.0) == A[1]
    @test etpg(10.1) == etpg(11) == etpg(148.298452) == A[end]

    etpf = @inferred(extrapolate(itpg, NaN))
    @test typeof(etpf) <: Interpolations.FilledExtrapolation
    @test parent(etpf) === itpg

    @test @inferred(size(etpf)) == (xmax,)
    @test isnan(@inferred(etpf(-2.5)))
    @test isnan(etpf(0.999))
    @test @inferred(etpf(1)) ≈ f(1)
    @test etpf(10) ≈ f(10)
    @test isnan(@inferred(etpf(10.001)))

    @test etpf(2.5,1) == etpf(2.5)   # for show method
    @test_throws BoundsError etpf(2.5,2)
    @test_throws BoundsError etpf(2.5,2,1)

    x =  @inferred(etpf(dual(-2.5,1)))
    @test isa(x, Dual)

    etpl = extrapolate(itpg, Line())
    k_lo = A[2] - A[1]
    x_lo = -3.2
    @test etpl(x_lo) ≈ A[1] + k_lo * (x_lo - 1)
    k_hi = A[end] - A[end-1]
    x_hi = xmax + 5.7
    @test etpl(x_hi) ≈ A[end] + k_hi * (x_hi - xmax)


    xmax, ymax = 8,8
    g(x, y) = (x^2 + 3x - 8) * (-2y^2 + y + 1)

    itp2g = interpolate(Float64[g(x,y) for x in 1:xmax, y in 1:ymax], (BSpline(Quadratic(Free(OnGrid()))), BSpline(Linear())))
    etp2g = extrapolate(itp2g, (Line(), Flat()))

    @test @inferred(etp2g(-0.5,4)) ≈ itp2g(1,4) - 1.5 * epsilon(etp2g(dual(1,1),4))
    @test @inferred(etp2g(5,100)) ≈ itp2g(5,ymax)

    etp2ud = extrapolate(itp2g, ((Line(), Flat()), Flat()))
    @test @inferred(etp2ud(-0.5,4)) ≈ itp2g(1,4) - 1.5 * epsilon(etp2g(dual(1,1),4))
    @test @inferred(etp2ud(5, -4)) == etp2ud(5,1)
    @test @inferred(etp2ud(100, 4)) == etp2ud(8,4)
    @test @inferred(etp2ud(-.5, 100)) == itp2g(1,8) - 1.5 * epsilon(etp2g(dual(1,1),8))

    etp2ll = extrapolate(itp2g, Line())
    @test @inferred(etp2ll(-0.5,100)) ≈ (itp2g(1,8) - 1.5 * epsilon(etp2ll(dual(1,1),8))) + (100 - 8) * epsilon(etp2ll(1,dual(8,1)))

    # Allow element types that don't support conversion to Int (#87):
    etp87g = extrapolate(interpolate([1.0im, 2.0im, 3.0im], BSpline(Linear())), 0.0im)
    @test @inferred(etp87g(1)) == 1.0im
    @test @inferred(etp87g(1.5)) == 1.5im
    @test @inferred(etp87g(0.75)) == 0.0im
    @test @inferred(etp87g(3.25)) == 0.0im

    # Make sure it works with Gridded too
    etp100g = extrapolate(interpolate(([10;20],),[100;110], Gridded(Linear())), Flat())
    @test @inferred(etp100g(5)) == 100
    @test @inferred(etp100g(15)) == 105
    @test @inferred(etp100g(25)) == 110
    # issue #178
    a = randn(10,10) + im*rand(10,10)
    etp = @inferred(extrapolate(interpolate((1:10, 1:10), a, Gridded(Linear())), 0.0))
    @test @inferred(etp(-1,0)) === 0.0+0.0im

    # check all extrapolations work with vectorized indexing
    for E in [0,Flat(),Line(),Periodic(),Reflect()]
        @test (@inferred(extrapolate(interpolate([0,0],BSpline(Linear())),E))([1.2, 1.8, 3.1])) == [0,0,0]
    end

    # Issue #156
    F     = *(collect(1.0:10.0), collect(1:4)')
    itp   = interpolate(F, (BSpline(Linear()), NoInterp()));
    itps   = scale(itp, 1:10, 1:4)
    itpe   = extrapolate(itps, (Line(), Throw()))
    @test itpe(10.1, 1) ≈ 10.1
    @test_throws BoundsError itpe(9.9, 0)

    # Issue #232
    targeterr = ArgumentError("cannot create a filled extrapolation with a type Line, use a value of this type instead (e.g., Line())")
    @test_throws targeterr extrapolate(itp, Line)

    # Issue #244
    xs = range(1e-2, stop = 8.3, length = 3)
    ys = sort(rand(3))
    itp = LinearInterpolation(xs, ys, extrapolation_bc = Flat())
    @test itp(8.3) ≈ ys[end]

    # Issue #288
    etp = extrapolate(interpolate(collect(1:12), NoInterp()) , 0)
    @test etp(0:23) == [0:12; zeros(Int, 11)]

    include("type-stability.jl")
    include("non1.jl")
end
