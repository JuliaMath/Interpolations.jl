module RunTests

using Base.Test
using Interpolations


# extrapolation tests
include("extrapolation/runtests.jl")
# b-spline interpolation tests
include("b-splines/runtests.jl")
include("nointerp.jl")

# scaling tests
include("scaling/runtests.jl")

# # test gradient evaluation
include("gradient.jl")

# gridded interpolation tests
include("gridded/runtests.jl")

# test interpolation with specific types
include("typing.jl")

# Tests copied from Grid.jl's old test suite
#include("grid.jl")

include("issues/runtests.jl")

include("io.jl")
include("convenience-constructors.jl")
include("readme-examples.jl")

end

nothing
