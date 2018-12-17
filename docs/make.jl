using Documenter, Interpolations
makedocs(
sitename="Interpolations.jl",
modules=[Interpolations],
format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
pages=["Home" => "index.md",
        "General usage" => "interpolations.md",
        "Interpolation algorithms" => "control.md",
        "Extrapolation" => "extrapolation.md",
        "Convenience Construcors" => "convenience-construction.md",
        "Library" => "api.md"]
)

deploydocs(repo="github.com/JuliaMath/Interpolations.jl")
