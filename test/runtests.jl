module RunTests

using Base.Test
using Interpolations

include("linear.jl")
include("quadratic.jl")

# Tests copied from Grid.jl's old test suite
include("grid.jl")

end
