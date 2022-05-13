using Unitful

@testset "Unitful scaled interpolation" begin
    # Abstract range
    Xs = (0:.2:5)u"m"
    Ys = Xs .^ 2
    sitp = scale(interpolate(Ys, BSpline(Quadratic(Free(OnGrid())))), Xs)
    @test sitp.((0:1:3)u"cm") ≈ ((0:1:3)u"cm") .^ 2

    # Step Range
    Xss = (0:5)u"m"
    Yss = Float64.(Xss .^ 2)
    sitp = scale(interpolate(Yss, BSpline(Quadratic(Free(OnGrid())))), Xss)
    @test sitp((0:1:3)u"cm") ≈ ((0:1:3)u"cm") .^ 2

end
