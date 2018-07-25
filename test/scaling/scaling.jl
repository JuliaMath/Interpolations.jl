module ScalingTests

using Interpolations, Compat
using Compat.Test, Compat.LinearAlgebra

# Model linear interpolation of y = -3 + .5x by interpolating y=x
# and then scaling to the new x range

itp = interpolate(1:1.0:10, BSpline(Linear()), OnGrid())

sitp = @inferred(scale(itp, -3:.5:1.5))
@test typeof(sitp) <: Interpolations.ScaledInterpolation
@test parent(sitp) === itp

for (x,y) in zip(-3:.05:1.5, 1:.1:10)
    @test sitp[x] ≈ y
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
    @test testfunction(x,y) ≈ sitp2[x,y]
end

# Test gradients of scaled grids
xs = -pi:.1:pi
ys = map(sin, xs)
itp = interpolate(ys, BSpline(Linear()), OnGrid())
sitp = @inferred scale(itp, xs)

for x in -pi:.1:pi
    g = @inferred(gradient(sitp, x))[1]
    @test ≈(cos(x),g,atol=0.05)
end

# Verify that return types are reasonable
@inferred(getindex(sitp2, -3.4, 1.2))
@inferred(getindex(sitp2, -3, 1))
@inferred(getindex(sitp2, -3.4, 1))

sitp32 = @inferred scale(interpolate(Float32[testfunction(x,y) for x in -5:.5:5, y in -4:.2:4], BSpline(Quadratic(Flat())), OnGrid()), -5f0:.5f0:5f0, -4f0:.2f0:4f0)
@test typeof(@inferred(getindex(sitp32, -3.4f0, 1.2f0))) == Float32

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
        dest[I] = sitp[I]
    end
    dest
end
rfoo = Array{Float64}(undef, Interpolations.ssize(sitp))
rbar = similar(rfoo)
foo!(rfoo, sitp)
bar!(rbar, sitp)
@test rfoo ≈ rbar

end
