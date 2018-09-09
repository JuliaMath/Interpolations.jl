using Interpolations, Test

@testset "MixedTests" begin
    A = rand(6,5)
    knots = (range(1, stop=size(A,1), length=size(A,1)), range(1, stop=size(A,2), length=size(A,2)))
    itp = @inferred(interpolate(knots, A, (Gridded(Linear()),NoInterp())))
    @inferred(itp(2, 2))
    @inferred(itp(CartesianIndex((2,2))))
    for j = 2:size(A,2)-1, i = 2:size(A,1)-1
        @test itp(i,j) ≈ A[i,j]
        @test itp(CartesianIndex(i,j)) ≈ A[i,j]
    end
    @inferred(itp(knots...)) ≈ A

    A = rand(8,20)
    knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
    itp = interpolate(knots, A, Gridded(Linear()))
    @test itp(4,1.2) ≈ A[2,6]

    @test_throws  ErrorException interpolate(knots, A, (Gridded(Linear()),NoInterp()))
end
