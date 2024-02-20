using Documenter
using FreeBird

makedocs(;
    pages=[
            "FreeBird.jl" => "index.md",
            "Potentials" => "Potentials.md",
            "EnergyEval" => "EnergyEval.md",
            "AbstractWalkers" => "AbstractWalkers.md",
            "AtomsMCMoves" => "AtomsMCMoves.md",
            "SamplingSchemes" => "SamplingSchemes.md",
            "API" => "API.md",
        ],
    sitename = "FreeBird.jl",
    format = Documenter.HTML(),
    modules = [FreeBird],
    doctest=false,
    checkdocs=:export,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/wexlergroup/FreeBird.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "feature/simple-lennard-jones", # just a test
)
