module ExtrapTests

using Base.Test
using Interpolations

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10
A = Float64[f(x) for x in 1:xmax]

itpg = interpolate(A, BSpline(Linear), OnGrid)

expg = extrapolate(itpg, Flat)

@test expg[-3] == expg[-4.5] == expg[0.9] == expg[1.0] == A[1]
@test expg[10.1] == expg[11] == expg[148.298452] == A[end]

end
