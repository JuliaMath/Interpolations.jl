if !isdefined(Main, :InterpolationTestUtils)
    include("InterpolationTestUtils.jl")
    @eval using Main.InterpolationTestUtils
end

using Test, SharedArrays, Random, ColorVectorSpace
using StaticArrays, WoodburyMatrices
using Interpolations


const isci = get(ENV, "CI", "") in ("true", "True")

@testset "Interpolations" begin
    @testset "method ambiguities" begin
        @test isempty(detect_ambiguities(Interpolations))
    end

    include("core.jl")
    # Hermite interpolation tests
    include("cubic_hermite.jl")
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

    # lanczos tests
    include("lanczos/runtests.jl")

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
    include("iterate.jl")
    include("pipes.jl")

    # Chain rules interaction
    include("chainrules.jl")

    include("gpu_support.jl")
end
