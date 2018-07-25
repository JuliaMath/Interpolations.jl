module ScalingWithExtrapTests

using Interpolations, Compat.Test
using Compat: range

xs = range(-5, stop=5, length=10)
ys = map(sin, xs)

function run_tests(sut::Interpolations.AbstractInterpolation{T,N,IT,OnGrid}, itp) where {T,N,IT}
    for x in xs
        @test ≈(sut[x],sin(x),atol=sqrt(eps(sin(x))))
    end
    @test sut[-5] == sut[-5.1] == sut[-15.8] == sut[-Inf] == itp[1]
    @test sut[5] == sut[5.1] == sut[15.8] == sut[Inf] == itp[end]
end

function run_tests(sut::Interpolations.AbstractInterpolation{T,N,IT,OnCell}, itp) where {T,N,IT}
    halfcell = (xs[2] - xs[1]) / 2

    for x in (5 + halfcell, 5 + 1.1halfcell, 15.8, Inf)
        @test sut[-x] == itp[.5]
        @test sut[x] == itp[end+.5]
    end
end

for GT in (OnGrid, OnCell)
    itp = interpolate(ys, BSpline(Quadratic(Flat())), GT())

    # Test extrapolating, then scaling
    eitp = extrapolate(itp, Flat())
    seitp = scale(eitp, xs)
    run_tests(seitp, itp)

    # Test scaling, then extrapolating
    sitp = scale(itp, xs)
    esitp = extrapolate(sitp, Flat())
    run_tests(esitp, itp)
end

end
