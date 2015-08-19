module ExtrapTests

using Base.Test
using Interpolations

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10
A = Float64[f(x) for x in 1:xmax]

itpg = interpolate(A, BSpline(Linear), OnGrid)

etpg = extrapolate(itpg, Flat)

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

etpf = @inferred(FilledInterpolation(itpg, 'x'))
@test_approx_eq etpf[2] f(2)
@test etpf[-1.5] == 'x'

end
