using Documenter, Interpolations
makedocs(
sitename="Interpolations.jl",
modules=[Interpolations],
pages=["Home" => "index.md",
        "General usage" => "interpolations.md",
        "Interpolation algorithms" => "control.md",
        "Extrapolation" => "extrapolation.md",
        "Convenience Construcors" => "convenience-construction.md",
        "Library" => "api.md"]
)
