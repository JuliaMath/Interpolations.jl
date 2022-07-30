using Documenter, Interpolations
makedocs(
sitename="Interpolations.jl",
modules=[Interpolations],
format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
pages=["Home" => "index.md",
        "Convenience Constructors" => "convenience-construction.md",
        "General usage" => "interpolations.md",
        "Interpolation algorithms" => "control.md",
        "Extrapolation" => "extrapolation.md",
        "Knot Iteration" => "iterate.md",
        "Developer documentation" => "devdocs.md",
        "Library" => "api.md",
        "News and Changes" => "NEWS.md",
        "Other Interpolation Packages" => "other_packages.md",
],
strict=true,
)

deploydocs(repo="github.com/JuliaMath/Interpolations.jl.git",
            push_preview=true)
