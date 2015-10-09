module TypingTests

using Interpolations, Base.Test

nx = 10
f(x) = convert(Float32, x^3/(nx-1))
g(x) = convert(Float32, 3x^2/(nx-1))

A = Float32[f(x) for x in 1:nx]

itp = interpolate(A, BSpline(Quadratic(Flat())), OnCell())

# display(plot(
#     layer(x=1:nx,y=[f(x) for x in 1:1//1:nx],Geom.point),
#     layer(x=1:.1:nx,y=[itp[x] for x in 1:1//10:nx],Geom.path),
# ))

for x in 3.1:.2:4.3
    @test_approx_eq_eps float(f(x)) float(itp[x]) abs(.1*f(x))
end

@test typeof(itp[3.5f0]) == Float32

for x in 3.1:.2:4.3
    @test_approx_eq_eps g(x) gradient(itp, x) abs(.1*g(x))
end

@test typeof(gradient(itp, 3.5f0)[1]) == Float32

# Rational element types
R = Rational{Int}[x^2//10 for x in 1:10]
itp = interpolate(R, BSpline(Quadratic(Free())), OnCell())
itp[11//10]

@test typeof(itp[11//10]) == Rational{Int}
@test itp[11//10] == (11//10)^2//10

end
