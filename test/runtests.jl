module RunTests

using Base.Test
using Interpolations

# b-spline interpolation tests
include("b-splines/runtests.jl")

# extrapolation tests
include("extrapolation/runtests.jl")

# # test gradient evaluation
# include("gradient.jl")

# # test interpolation with specific types
# include("typing.jl")

# Tests copied from Grid.jl's old test suite
#include("grid.jl")

println("all tests passed")

end
