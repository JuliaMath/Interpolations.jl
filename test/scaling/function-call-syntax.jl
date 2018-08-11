module ScalingFunctionCallTests

using Interpolations, Test

# Model linear interpolation of y = -3 + .5x by interpolating y=x
# and then scaling to the new x range

itp = interpolate(1:1.0:10, BSpline(Linear()), OnGrid())

sitp = @inferred(scale(itp, -3:.5:1.5))
@test typeof(sitp) <: Interpolations.ScaledInterpolation
@test parent(sitp) === itp

for (x,y) in zip(-3:.05:1.5, 1:.1:10)
    @test sitp(x) ≈ y
end

# Verify that it works in >1D, with different types of ranges

gauss(phi, mu, sigma) = exp(-(phi-mu)^2 / (2sigma)^2)
testfunction(x,y) = gauss(x, 0.5, 4) * gauss(y, -.5, 2)

xs = -5:.5:5
ys = -4:.2:4
zs = Float64[testfunction(x,y) for x in xs, y in ys]

itp2 = interpolate(zs, BSpline(Quadratic(Flat())), OnGrid())
sitp2 = @inferred scale(itp2, xs, ys)

for x in xs, y in ys
    @test testfunction(x,y) ≈ sitp2(x,y)
end

# Iteration
itp = interpolate(rand(3,3,3), BSpline(Quadratic(Flat())), OnCell())
knots = map(d->1:10:21, 1:3)
sitp = @inferred scale(itp, knots...)

iter = @inferred(eachvalue(sitp))

@static if VERSION < v"0.7.0-DEV.5126"
    state = @inferred(start(iter))
    @test !(@inferred(done(iter, state)))
    val, state = @inferred(next(iter, state))
else
    iter_next = iterate(iter)
    @test iter_next isa Tuple
    @test iter_next[1] isa Float64
    state = iter_next[2]
    inferred_next = Base.return_types(iterate, (typeof(iter),))
    @test length(inferred_next) == 1
    @test inferred_next[1] == Union{Nothing,Tuple{Float64,typeof(state)}}
    iter_next = iterate(iter, state)
    @test iter_next isa Tuple
    @test iter_next[1] isa Float64
    inferred_next = Base.return_types(iterate, (typeof(iter),typeof(state)))
    state = iter_next[2]
    @test length(inferred_next) == 1
    @test inferred_next[1] == Union{Nothing,Tuple{Float64,typeof(state)}}
end

function foo!(dest, sitp)
    i = 0
    for s in eachvalue(sitp)
        dest[i+=1] = s
    end
    dest
end
function bar!(dest, sitp)
    for I in CartesianIndices(size(dest))
        dest[I] = sitp(I)
    end
    dest
end
rfoo = Array{Float64}(undef, Interpolations.ssize(sitp))
rbar = similar(rfoo)
foo!(rfoo, sitp)
bar!(rbar, sitp)
@test rfoo ≈ rbar

# with extrapolation
END = 10
xs = range(-5, stop=5, length=END)
ys = map(sin, xs)

function run_tests(sut::Interpolations.AbstractInterpolation{T,N,IT,OnGrid}, itp) where {T,N,IT}
    for x in xs
        @test ≈(sut[x],sin(x),atol=sqrt(eps(sin(x))))
    end
    @test sut(-5) == sut(-5.1) == sut(-15.8) == sut(-Inf) == itp(1)
    @test sut(5) == sut(5.1) == sut(15.8) == sut(Inf) == itp(END)
end

function run_tests(sut::Interpolations.AbstractInterpolation{T,N,IT,OnCell}, itp) where {T,N,IT}
    halfcell = (xs[2] - xs[1]) / 2

    for x in (5 + halfcell, 5 + 1.1halfcell, 15.8, Inf)
        @test sut(-x) == itp(.5)
        @test sut(x) == itp(END+.5)
    end
end

for GT in (OnGrid, OnCell)
    itp3 = interpolate(ys, BSpline(Quadratic(Flat())), GT())

    # Test extrapolating, then scaling
    eitp = extrapolate(itp3, Flat())
    seitp = scale(eitp, xs)
    run_tests(seitp, itp3)

    # Test scaling, then extrapolating
    sitp3 = scale(itp3, xs)
    esitp = extrapolate(sitp3, Flat())
    run_tests(esitp, itp3)
end

end
