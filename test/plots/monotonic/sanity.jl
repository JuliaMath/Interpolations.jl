using Plots
using Interpolations

pyplot()

x = [0.0, 0.2, 0.5, 0.6, 0.9, 1.0]
y = [10.0, 0.0, -5.0, 10.0, -8.0, -2.0]
grid = 0.0:0.002:1.0

itypes = [LinearMonotonicInterpolation(),
    FiniteDifferenceMonotonicInterpolation(),
    CardinalMonotonicInterpolation(0.0),
    CardinalMonotonicInterpolation(0.5),
    CardinalMonotonicInterpolation(1.0),
    FritschCarlsonMonotonicInterpolation(),
    FritschButlandMonotonicInterpolation(),
    SteffenMonotonicInterpolation()
    ]

scatter(x, y, label = "data")
for itype in itypes
    itp = interpolate(x, y, itype)
    plot!([x for x in grid], [itp(x) for x in grid], label = string(itype))
end
gui()
