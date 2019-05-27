@testset "BSpline" begin
    include("constant.jl")
    include("linear.jl")
    include("quadratic.jl")
    include("cubic.jl")
    include("mixed.jl")
    include("multivalued.jl")
    include("non1.jl")

    @test eltype(@inferred(interpolate(rand(Float16, 3, 3), BSpline(Linear())))) == Float16 # issue #308
end
