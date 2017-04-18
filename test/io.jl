module IOTests

using Base.Test
using Interpolations

SPACE = if VERSION < v"0.6.0-dev.2505" # julia PR #20288
    ""
else
    " "
end

@testset "BSpline" begin
    A = rand(8,20)

    itp = interpolate(A, BSpline(Constant()), OnCell())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Constant()), OnCell()) with element type Float64"

    itp = interpolate(A, BSpline(Constant()), OnGrid())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Constant()), OnGrid()) with element type Float64"

    itp = interpolate(A, BSpline(Linear()), OnGrid())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Linear()), OnGrid()) with element type Float64"

    itp = interpolate(A, BSpline(Quadratic(Reflect())), OnCell())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Quadratic(Reflect())), OnCell()) with element type Float64"

    itp = interpolate(A, (BSpline(Linear()), NoInterp()), OnGrid())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, (BSpline(Linear()), NoInterp()), OnGrid()) with element type Float64"

    itp = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    @test summary(itp) == "8×20 interpolate(::Array{Float64,2}, BSpline(Quadratic(InPlace())), OnCell()) with element type Float64"
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
end

@testset "scaled" begin
    itp = interpolate(1:1.0:10, BSpline(Linear()), OnGrid())
    sitp = scale(itp, -3:.5:1.5)
    @test summary(sitp) == "10-element scale(interpolate(::Array{Float64,1}, BSpline(Linear()), OnGrid()), (-3.0:0.5:1.5,)) with element type Float64"

    gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
    testfunction(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)
    xs = -5:.5:5
    ys = -4:.2:4
    zs = Float64[testfunction(x,y) for x in xs, y in ys]
    itp2 = interpolate(zs, BSpline(Quadratic(Flat())), OnGrid())
    sitp2 = scale(itp2, xs, ys)
    @test summary(sitp2) == "21×41 scale(interpolate(::Array{Float64,2}, BSpline(Quadratic(Flat())), OnGrid()), (-5.0:0.5:5.0,$SPACE-4.0:0.2:4.0)) with element type Float64"
end

@testset "Extrapolation" begin
    A = rand(8,20)

    itpg = interpolate(A, BSpline(Linear()), OnGrid())
    etpg = extrapolate(itpg, Flat())
    @test summary(etpg) == "8×20 extrapolate(interpolate(::Array{Float64,2}, BSpline(Linear()), OnGrid()), Flat()) with element type Float64"

    etpf = extrapolate(itpg, NaN)
    @test summary(etpf) == "8×20 extrapolate(interpolate(::Array{Float64,2}, BSpline(Linear()), OnGrid()), NaN) with element type Float64"
end


end # Module
