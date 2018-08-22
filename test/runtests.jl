if !isdefined(Main, :InterpolationTestUtils)
    include("InterpolationTestUtils.jl")
    @eval using Main.InterpolationTestUtils
end

using Test, SharedArrays, Random
using StaticArrays, WoodburyMatrices
ambs = detect_ambiguities(StaticArrays, WoodburyMatrices, Base, Core)

using Interpolations
@test isempty(setdiff(detect_ambiguities(Interpolations, Base, Core), ambs))

@testset "Interpolations" begin

    # b-spline interpolation tests
    include("b-splines/runtests.jl")
    include("nointerp.jl")
    # extrapolation tests
    include("extrapolation/runtests.jl")

    # scaling tests
    include("scaling/runtests.jl")

    # # test gradient evaluation
    include("gradient.jl")
    # # test hessian evaluation
    include("hessian.jl")

    # # gridded interpolation tests
    # include("gridded/runtests.jl")

    # test interpolation with specific types
    include("typing.jl")

    # Tests copied from Grid.jl's old test suite
    # include("grid.jl")

    include("issues/runtests.jl")

    include("io.jl")
    # include("convenience-constructors.jl")
    include("readme-examples.jl")

end
