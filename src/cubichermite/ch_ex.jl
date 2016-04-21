using Interpolations
include("cubichermite.jl")


# Use points 1, 2, 3, 4, 5
x = collect(1:5)

# Test with exponential function
f(x) = log(x)
fp(x) = 1.0 ./ x

# Create arrays we need
y = f(x)
m = 1.0 ./ fp(x)

# Create the cubic hermite interpolation coeffs
# Note: m is the derivatives, so this has more info
c = pre_solve(y, m)
c_nd = pre_solve(y)

# Create an interpolations cubic b-spline
itp = interpolate(y, BSpline(Cubic(Natural())), OnGrid())

# Create arrays to fill
N = 250
out_chp = Array(Float64, N)
out_chp_nd = Array(Float64, N)
out_itp = Array(Float64, N)
out_tru = Array(Float64, N)

for (i, v) in enumerate(linspace(1.01, 4.99, N))
    out_chp[i] = evaluate(c, v)
    out_chp_nd[i] = evaluate(c_nd, v)
    out_itp[i] = itp[v]
    out_tru[i] = f(v)
end

println("Max distance from cubic hermite (derivative) to truth is ", maxabs(out_chp-out_tru))
println("Max distance from cubic hermite (no derivative) to truth is ", maxabs(out_chp_nd-out_tru))
println("Max distance from cubic Bspline to truth is ", maxabs(out_itp-out_tru))

