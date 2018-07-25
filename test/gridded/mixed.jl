module MixedTests

using Interpolations, Compat.Test
using Compat: range

A = rand(6,5)
knots = (collect(range(1, stop=size(A,1), length=size(A,1))),collect(range(1, stop=size(A,2), length=size(A,2))))
itp = @inferred(interpolate(knots, A, (Gridded(Linear()),NoInterp())))
@inferred(getindex(itp, 2, 2))
@inferred(getindex(itp, CartesianIndex((2,2))))
for j = 2:size(A,2)-1, i = 2:size(A,1)-1
    @test itp[i,j] ≈ A[i,j]
    @test itp[CartesianIndex((i,j))] ≈ A[i,j]
end
@test itp[knots...] ≈ A
@inferred(getindex(itp, knots...))

A = rand(8,20)
knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
itp = interpolate(knots, A, Gridded(Linear()))
@test itp[4,1.2] ≈ A[2,6]

@test_throws  ErrorException interpolate(knots, A, (Gridded(Linear()),NoInterp()))
end
