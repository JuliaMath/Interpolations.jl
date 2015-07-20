
module ScalingWithExtrapTests

using Interpolations, Base.Test

xs = linspace(-5, 5, 10)
ys = sin(xs)

function run_tests{T,N,IT}(sut::Interpolations.AbstractInterpolation{T,N,IT,OnGrid}, itp)
    for x in xs
        @test_approx_eq_eps sut[x] sin(x) sqrt(eps(sin(x)))
    end
    @test sut[-5] == sut[-5.1] == sut[-15.8] == sut[-Inf] == itp[1]
    @test sut[5] == sut[5.1] == sut[15.8] == sut[Inf] == itp[end]
end

function run_tests{T,N,IT}(sut::Interpolations.AbstractInterpolation{T,N,IT,OnCell}, itp)
    halfcell = (xs[2] - xs[1]) / 2

    for x in (5 + halfcell, 5 + 1.1halfcell, 15.8, Inf)
        @test sut[-x] == itp[.5]
        @test sut[x] == itp[end+.5]
    end
end

for GT in (OnGrid, OnCell)
    itp = interpolate(ys, BSpline(Quadratic(Flat)), GT)

    # Test extrapolating, then scaling
    eitp = extrapolate(itp, Flat)
    seitp = scale(eitp, xs)
    run_tests(seitp, itp)

    # Test scaling, then extrapolating
    sitp = scale(itp, xs)
    esitp = extrapolate(sitp, Flat)
    run_tests(esitp, itp)
end

end