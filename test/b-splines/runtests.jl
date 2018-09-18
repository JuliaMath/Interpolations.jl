@testset "BSpline" begin
    include("constant.jl")
    include("linear.jl")
    include("quadratic.jl")
    include("cubic.jl")
    include("mixed.jl")
    include("multivalued.jl")
    include("non1.jl")
end
