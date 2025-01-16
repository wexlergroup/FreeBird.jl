using Documenter
using FreeBird
using Plots

makedocs(;
    pages=[
            "FreeBird.jl" => "index.md",
            "Quick Start" => "quick_start.md",
            "Modules and Documentation" => [
            "AbstractLiveSets" => "AbstractLiveSets.md",
            "AbstractWalkers" => "AbstractWalkers.md",
            "SamplingSchemes" => "SamplingSchemes.md",
            "AbstractPotentials" => "AbstractPotentials.md",
            "AbstractHamiltonians" => "AbstractHamiltonians.md",
            "EnergyEval" => "EnergyEval.md",
            "MonteCarloMoves" => "MonteCarloMoves.md",
            "FreeBirdIO" => "FreeBirdIO.md",
            "AnalysisTools" => "AnalysisTools.md"],
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
    devbranch = "dev",
    push_preview = true,
)
