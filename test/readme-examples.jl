# verify examples from README.md run
using Interpolations, Test

@testset "Readme Examples" begin

    ## Bsplines
    a = randn(6)
    A = randn(5, 5)

    # Nearest-neighbor interpolation
    itp = interpolate(a, BSpline(Constant()))
    v = itp(5.4)   # returns a[5]
    @test v ≈ a[5]

    # (Multi)linear interpolation
    itp = interpolate(A, BSpline(Linear()))
    v = itp(3.2, 4.1)  # returns 0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5])
    @test v ≈ (0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5]))

    # Quadratic interpolation with reflecting boundary conditions
    # Quadratic is the lowest order that has continuous gradient
    itp = interpolate(A, BSpline(Quadratic(Reflect(OnCell()))))

    # Linear interpolation in the first dimension, and no interpolation (just lookup) in the second
    itp = interpolate(A, (BSpline(Linear()), NoInterp()))
    v = itp(3.65, 5)  # returns  0.35*A[3,5] + 0.65*A[4,5]
    @test v ≈ (0.35*A[3,5] + 0.65*A[4,5])


    ## Scaled Bsplines
    A_x = 1.:2.:40.
    A = [log(x) for x in A_x]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, A_x)
    @test sitp(3.) ≈ log(3.) # exactly log(3.)
    @test sitp(3.5) ≈ log(3.5) atol=.1 # approximately log(3.5)

    # For multidimensional uniformly spaced grids
    A_x1 = 1:.1:10
    A_x2 = 1:.5:20
    f(x1, x2) = log(x1+x2)
    A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, A_x1, A_x2)
    @test sitp(5., 10.) ≈ log(5 + 10) # exactly log(5 + 10)
    @test sitp(5.6, 7.1) ≈ log(5.6 + 7.1) atol=.1 # approximately log(5.6 + 7.1)

    ## Gridded interpolation
    A = rand(8,20)
    knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
    itp = interpolate(knots, A, Gridded(Linear()))
    @test itp(4,1.2) ≈ A[2,6] atol=.1  # approximately A[2,6]

end
