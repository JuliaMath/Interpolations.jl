module GradientTests
println("Testing gradient evaluation")
using Base.Test, Interpolations

nx = 10
f1(x) = sin((x-3)*2pi/(nx-1) - 1)
g1(x) = 2pi/(nx-1) * cos((x-3)*2pi/(nx-1) - 1)

# Gradient of Constant should always be 0
itp1 = Interpolation(Float64[f1(x) for x in 1:nx-1],
            Constant(OnGrid()), ExtrapPeriodic())
for x in 1:nx
    @test gradient(itp1, x)[1] == 0
end

# Since Linear is OnGrid in the domain, check the gradients between grid points
itp1 = Interpolation(Float64[f1(x) for x in 1:nx-1],
            Linear(OnGrid()), ExtrapPeriodic())
for x in 2.5:nx-1.5
    @test_approx_eq_eps g1(x) gradient(itp1, x)[1] abs(.1*g1(x))
end

# Since Quadratic is OnCell in the domain, check gradients at grid points
itp1 = Interpolation(Float64[f1(x) for x in 1:nx-1], 
            Quadratic(Periodic(),OnGrid()), ExtrapPeriodic())
for x in 2:nx-1
    @test_approx_eq_eps g1(x) gradient(itp1, x)[1] abs(.05*g1(x))
end

end
