using Documenter
using FreeBird

makedocs(
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
)
