using Zygote

x = collect(1:10)
y = sin.(x)

# Simple test that rrules give equivalent results to internal gradient function

itp = interpolate(y,BSpline(Linear()))

@test Zygote.gradient(itp, 1)[1] == Interpolations.gradient(itp, 1)