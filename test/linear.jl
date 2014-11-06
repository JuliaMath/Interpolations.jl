
## Tests for linear interpolation of scalars

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

A = Float64[f(x) for x in 1:xmax]

## ExtrapError

itp1 = Interpolation(A, Linear, ExtrapError)

for x in [3.1:.2:4.3]  
    @test_approx_eq_eps f(x) itp1[x] abs(.1*f(x))
end

@test_throws BoundsError itp1[-3]