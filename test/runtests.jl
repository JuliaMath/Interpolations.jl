if !isdefined(Main, :InterpolationTestUtils)
    include("InterpolationTestUtils.jl")
    @eval using Main.InterpolationTestUtils
end

using Test, SharedArrays, Random
using StaticArrays, WoodburyMatrices
ambs = detect_ambiguities(StaticArrays, WoodburyMatrices, Base, Core)

using Interpolations
@test isempty(setdiff(detect_ambiguities(Interpolations, Base, Core), ambs))

const isci = get(ENV, "CI", "") in ("true", "True")

@testset "Interpolations" begin
    include("core.jl")

    # b-spline interpolation tests
    include("b-splines/runtests.jl")
    isci && println("finished b-spline")
    include("nointerp.jl")
    # extrapolation tests
    include("extrapolation/runtests.jl")

    # scaling tests
    include("scaling/runtests.jl")

    # monotonic tests
    include("monotonic/runtests.jl")

    # test gradient evaluation
    include("gradient.jl")
    isci && println("finished gradient")
    # test hessian evaluation
    include("hessian.jl")
    isci && println("finished hessian")

    # gridded interpolation tests
    include("gridded/runtests.jl")

    # test interpolation with specific types
    include("typing.jl")

    include("issues/runtests.jl")

    include("io.jl")
    include("convenience-constructors.jl")
    include("readme-examples.jl")

end
