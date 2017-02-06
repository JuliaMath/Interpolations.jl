module ExtrapTests

using Base.Test, DualNumbers
using Interpolations

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10
A = Float64[f(x) for x in 1:xmax]

itpg = interpolate(A, BSpline(Linear()), OnGrid())

etpg = extrapolate(itpg, Flat())

@test etpg[-3] == etpg[-4.5] == etpg[0.9] == etpg[1.0] == A[1]
@test etpg[10.1] == etpg[11] == etpg[148.298452] == A[end]

etpf = @inferred(extrapolate(itpg, NaN))

@test @inferred(size(etpf)) == (xmax,)
@test isnan(@inferred(getindex(etpf, -2.5)))
@test isnan(etpf[0.999])
@test_approx_eq @inferred(getindex(etpf, 1)) f(1)
@test_approx_eq etpf[10] f(10)
@test isnan(@inferred(getindex(etpf,10.001)))

@test etpf[2.5,1] == etpf[2.5]   # for show method
@test_throws BoundsError etpf[2.5,2]
@test_throws ErrorException etpf[2.5,2,1]  # this will probably become a BoundsError someday

@test isa(@inferred(getindex(etpf, dual(-2.5,1))), Dual)

etpl = extrapolate(itpg, Linear())
k_lo = A[2] - A[1]
x_lo = -3.2
@test_approx_eq etpl[x_lo] A[1]+k_lo*(x_lo-1)
k_hi = A[end] - A[end-1]
x_hi = xmax + 5.7
@test_approx_eq etpl[x_hi] A[end]+k_hi*(x_hi-xmax)


xmax, ymax = 8,8
g(x, y) = (x^2 + 3x - 8) * (-2y^2 + y + 1)

itp2g = interpolate(Float64[g(x,y) for x in 1:xmax, y in 1:ymax], (BSpline(Quadratic(Free())), BSpline(Linear())), OnGrid())
etp2g = extrapolate(itp2g, (Linear(), Flat()))

@test_approx_eq @inferred(getindex(etp2g, -.5, 4)) itp2g[1,4]-1.5*epsilon(etp2g[dual(1,1),4])
@test_approx_eq @inferred(getindex(etp2g, 5, 100)) itp2g[5,ymax]

etp2ud = extrapolate(itp2g, ((Linear(), Flat()), Flat()))
@test_approx_eq @inferred(getindex(etp2ud, -.5, 4)) itp2g[1,4] - 1.5*epsilon(etp2g[dual(1,1),4])
@test @inferred(getindex(etp2ud, 5, -4)) == etp2ud[5,1]
@test @inferred(getindex(etp2ud, 100, 4)) == etp2ud[8,4]
@test @inferred(getindex(etp2ud, -.5, 100)) == itp2g[1,8] - 1.5 * epsilon(etp2g[dual(1,1),8])

etp2ll = extrapolate(itp2g, Linear())
@test_approx_eq @inferred(getindex(etp2ll, -.5, 100)) itp2g[1,8]-1.5*epsilon(etp2ll[dual(1,1),8]) + (100-8) * epsilon(etp2ll[1,dual(8,1)])

# Allow element types that don't support conversion to Int (#87):
etp87g = extrapolate(interpolate([1.0im, 2.0im, 3.0im], BSpline(Linear()), OnGrid()), 0.0im)
@test @inferred(getindex(etp87g, 1)) == 1.0im
@test @inferred(getindex(etp87g, 1.5)) == 1.5im
@test @inferred(getindex(etp87g, 0.75)) == 0.0im
@test @inferred(getindex(etp87g, 3.25)) == 0.0im

etp87c = extrapolate(interpolate([1.0im, 2.0im, 3.0im], BSpline(Linear()), OnCell()), 0.0im)
@test @inferred(getindex(etp87c, 1)) == 1.0im
@test @inferred(getindex(etp87c, 1.5)) == 1.5im
@test @inferred(getindex(etp87c, 0.75)) == 0.75im
@test @inferred(getindex(etp87c, 3.25)) == 3.25im
@test @inferred(getindex(etp87g, 0)) == 0.0im
@test @inferred(getindex(etp87g, 3.7)) == 0.0im

# Make sure it works with Gridded too
etp100g = extrapolate(interpolate(([10;20],),[100;110], Gridded(Linear())), Flat())
@test @inferred(getindex(etp100g, 5)) == 100
@test @inferred(getindex(etp100g, 15)) == 105
@test @inferred(getindex(etp100g, 25)) == 110
end

include("type-stability.jl")