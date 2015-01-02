module RunTests

using Base.Test
using Interpolations

# specific tests for linear extrapolation
include("linear.jl")
# specific tests for quadratic extrapolation
include("quadratic.jl")

# tests that ensure that itp[i...] = A[i...] for integer
# indices inbounds in A.
include("on-grid.jl")

# test gradient evaluation
include("gradient.jl")

# Tests copied from Grid.jl's old test suite
#include("grid.jl")

end
