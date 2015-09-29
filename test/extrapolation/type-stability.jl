module ExtrapTypeStability

using Base.Test, Interpolations, DualNumbers

# Test type-stability of 1-dimensional extrapolation
f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10
A = Float64[f(x) for x in 1:xmax]
itpg = interpolate(A, BSpline(Linear), OnGrid)

schemes = (
    Flat,
    Linear,
    Reflect,
    Periodic
)

for etp in map(E -> extrapolate(itpg, E), schemes),
    x in [
        # In-bounds evaluation
        3.4, 3, dual(3.1),
        # Out-of-bounds evaluation
        -3.4, -3, dual(-3,1),
        13.4, 13, dual(13,1)
    ]
    @inferred(getindex(etp, x))
end

# Test type-stability of 2-dimensional extrapolation with homogeneous scheme
g(y) = (y/100)^3
ymax = 4
A = Float64[f(x)*g(y) for x in 1:xmax, y in 1:ymax]
itp2 = interpolate(A, BSpline(Linear), OnGrid)

for (etp2,E) in map(E -> (extrapolate(itp2, E), E), schemes),
    x in (
        # In-bounds evaluation
        3.4, 3, dual(3.1),
        # Out-of-bounds evaluation
        -3.4, -3, dual(-3,1),
        13.4, 13, dual(13,1)
    ),
    y in (
        # In-bounds evaluation
        2.1, 2, dual(2.3, 1),
        # Out-of-bounds evaluation
        -2.1, -2, dual(-2.3, 1),
        12.1, 12, dual(12.1, 1)
    )
    @inferred(getindex(etp2, x, y))
end

end