module TypingTests

using Interpolations, Base.Test

nx = 10
f(x) = convert(Float32, x^3/(nx-1))
g(x) = convert(Float32, 3x^2/(nx-1))

A = Float32[f(x) for x in 1:nx]

itp = Interpolation(A, Quadratic(Flat(), OnCell()), ExtrapConstant())

# display(plot(
#     layer(x=1:nx,y=[f(x) for x in 1:1//1:nx],Geom.point),
#     layer(x=1:.1:nx,y=[itp[x] for x in 1:1//10:nx],Geom.path),
# ))

for x in [3.1:.2:4.3]
    @test_approx_eq_eps float(f(x)) float(itp[x]) abs(.1*f(x))
end

@test typeof(itp[3.5]) == Float32

for x in [3.1:.2:4.3]
    @test_approx_eq_eps g(x) gradient(itp, x) abs(.1*g(x))
end

@test typeof(gradient(itp, 3.5)[1]) == Float32

end
