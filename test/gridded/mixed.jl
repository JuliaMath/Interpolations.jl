module MixedTests

using Interpolations, Base.Test

A = rand(6,5)
knots = (collect(linspace(1,size(A,1),size(A,1))),collect(linspace(1,size(A,2),size(A,2))))
itp = @inferred(interpolate(knots, A, (Gridded(Linear()),NoInterp())))
@inferred(getindex(itp, 2, 2))
@inferred(getindex(itp, CartesianIndex((2,2))))
for j = 2:size(A,2)-1, i = 2:size(A,1)-1
    @test_approx_eq itp[i,j] A[i,j]
    @test_approx_eq itp[CartesianIndex((i,j))] A[i,j]
end
@test_approx_eq itp[knots...] A
@inferred(getindex(itp, knots...))

A = rand(8,20)
knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
itp = interpolate(knots, A, Gridded(Linear()))
@test_approx_eq itp[4,1.2] A[2,6]

@test_throws  ErrorException interpolate(knots, A, (Gridded(Linear()),NoInterp()))
end
