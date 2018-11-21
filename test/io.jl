using Interpolations
using Test

@testset "IO" begin
    SPACE = " "

    @testset "BSpline" begin
        A = rand(8,20)

        itp = interpolate(A, BSpline(Constant()))
        @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Constant())) with element type Float64"

        itp = interpolate(A, BSpline(Constant()))
        @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Constant())) with element type Float64"

        itp = interpolate(A, BSpline(Linear()))
        @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Linear())) with element type Float64"

        itp = interpolate(A, BSpline(Quadratic(Reflect(OnCell()))))
        @test summary(itp) == "8×20 interpolate(OffsetArray(::Array{Float64,2}, 0:9, 0:21), BSpline(Quadratic(Reflect(OnCell())))) with element type Float64"

        itp = interpolate(A, (BSpline(Linear()), NoInterp()))
        @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, (BSpline(Linear()), NoInterp())) with element type Float64"

        itp = interpolate!(copy(A), BSpline(Quadratic(InPlace(OnCell()))))
        @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Quadratic(InPlace(OnCell())))) with element type Float64"
    end

    @testset "Gridded" begin
        A = rand(20)
        A_x = collect(1.0:2.0:40.0)
        knots = (A_x,)
        itp = interpolate(knots, A, Gridded(Linear()))
        @test summary(itp) == "20-element interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Linear())) with element type Float64"

        A = rand(8,20)
        knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
        itp = interpolate(knots, A, Gridded(Linear()))
        @test summary(itp) == "8×20 interpolate((::Array{Int64,1},::Array{Float64,1}), ::Array{Float64,2}, Gridded(Linear())) with element type Float64"

        itp = interpolate(knots, A, (Gridded(Linear()),Gridded(Constant())))
        @test summary(itp) == "8×20 interpolate((::Array{Int64,1},::Array{Float64,1}), ::Array{Float64,2}, (Gridded(Linear()), Gridded(Constant()))) with element type Float64"

        # issue #260
        A = (1:4)/4
        itp = interpolate((range(0.0, stop=0.3, length=4),), A, Gridded(Linear()))
        io = IOBuffer()
        show(io, MIME("text/plain"), itp)
        str = String(take!(io))
        @test str == "4-element interpolate((0.0:0.1:0.3,), ::Array{Float64,1}, Gridded(Linear())) with element type Float64:\n 0.25\n 0.5 \n 0.75\n 1.0 "
    end

    @testset "scaled" begin
        itp = interpolate(1:1.0:10, BSpline(Linear()))
        sitp = scale(itp, -3:.5:1.5)
        @test summary(sitp) == "10-element scale(interpolate(::Array{Float64,1}, BSpline(Linear())), (-3.0:0.5:1.5,)) with element type Float64"
        io = IOBuffer()
        show(io, MIME("text/plain"), sitp)
        str = String(take!(io))
        @test str == "10-element scale(interpolate(::Array{Float64,1}, BSpline(Linear())), (-3.0:0.5:1.5,)) with element type Float64:\n  1.0\n  2.0\n  3.0\n  4.0\n  5.0\n  6.0\n  7.0\n  8.0\n  9.0\n 10.0"

        gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
        testfunction(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)
        xs = -5:.5:5
        ys = -4:.2:4
        zs = Float64[testfunction(x,y) for x in xs, y in ys]
        itp2 = interpolate(zs, BSpline(Quadratic(Flat(OnGrid()))))
        sitp2 = scale(itp2, xs, ys)
        @test summary(sitp2) == "21×41 scale(interpolate(OffsetArray(::Array{Float64,2}, 0:22, 0:42), BSpline(Quadratic(Flat(OnGrid())))), (-5.0:0.5:5.0,$SPACE-4.0:0.2:4.0)) with element type Float64"
    end

    @testset "Monotonic" begin
        A = rand(20)
        A_x = collect(1.0:2.0:40.0)
        itp = interpolate(A_x, A, FritschButlandMonotonicInterpolation())
        @test summary(itp) == "20-element interpolate(::Array{Float64,1}, ::Array{Float64,1}, FritschButlandMonotonicInterpolation()) with element type Float64"
    end

    @testset "Extrapolation" begin
        A = rand(8,20)

        itpg = interpolate(A, BSpline(Linear()))
        etpg = extrapolate(itpg, Flat())
        @test summary(etpg) == "8×20 extrapolate(interpolate(::Array{Float64,2}, BSpline(Linear())), Flat()) with element type Float64"

        etpf = extrapolate(itpg, NaN)
        @test summary(etpf) == "8×20 extrapolate(interpolate(::Array{Float64,2}, BSpline(Linear())), NaN) with element type Float64"
    end

    @testset "Combinations" begin
        xs = 1.0:0.5:2.0
        ys = xs.^3
        interp_cubic = CubicSplineInterpolation(xs, ys)
        io = IOBuffer()
        show(io, MIME("text/plain"), interp_cubic)
        str = String(take!(io))
        @test occursin("3.375", str)
    end

    @testset "Issue #260" begin
        itp = scale(interpolate(zeros(3,3,3),BSpline(Cubic(Free(OnGrid())))),2:4,2:4,2:4)
        io = IOBuffer()
        show(io, MIME("text/plain"), itp)
        str = String(take!(io))
        @test str == "3×3×3 scale(interpolate(OffsetArray(::Array{Float64,3}, 0:4, 0:4, 0:4), BSpline(Cubic(Free(OnGrid())))), (2:4, 2:4, 2:4)) with element type Float64:\n[:, :, 1] =\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n\n[:, :, 2] =\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n\n[:, :, 3] =\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n 0.0  0.0  0.0"
    end
end # testset
