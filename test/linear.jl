module Linear1DTests
println("Testing Linear interpolation in 1D...")
using Interpolations, Test

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

A = Float64[f(x) for x in 1:xmax]

## ExtrapError

itp1 = Interpolation(A, Linear(OnGrid()), ExtrapError())

for x in 1:.2:xmax
    @test ≈(f(x),itp1[x],atol=abs(0.1 * f(x)))
end

@test_throws BoundsError itp1[-3]
@test_throws BoundsError itp1[xmax+.1]

## ExtrapNaN

itp2 = Interpolation(A, Linear(OnGrid()), ExtrapNaN())

for x in 1:.2:xmax
    @test ≈(f(x),itp2[x],atol=abs(0.1 * f(x)))
end

xlo, xhi = itp2[.9], itp2[xmax+.2]
@test isnan(xlo)
@test isnan(xhi)

## ExtrapConstant

itp3 = Interpolation(A, Linear(OnGrid()), ExtrapConstant())
for x in 3.1:.2:4.3
    @test ≈(f(x),itp3[x],atol=abs(0.1 * f(x)))
end

xlo, xhi = itp3[.9], itp3[xmax+.2]
@test xlo == A[1]
@test xhi == A[end]

# Rational element types
A = Rational{Int}[x//10+2 for x in 1:10]
itp = Interpolation(A, Linear(OnGrid()), ExtrapNaN())
@test itp[11//10] == (11//10)//10+2

end

module Linear2DTests
println("Testing Linear interpolation in 2D...")
using Interpolations, Base.Test

f(x,y) = sin(x/10)*cos(y/6)

xmax, ymax = 30,10

A = Float64[f(x,y) for x in 1:xmax, y in 1:ymax]

itp1 = Interpolation(A, Linear(OnGrid()), ExtrapError())

for x in 2.1:.2:xmax-1, y in 1.9:.2:ymax-.9
    @test ≈(f(x,y),itp1[x,y],atol=0.01)
end

end
