using Documenter
using Test

@testset "Doctests" begin
    DocMeta.setdocmeta!(
        Interpolations,
        :DocTestSetup,
        :(using Interpolations),
        recursive=true,
    )
    doctest(Interpolations, manual=true)
end
