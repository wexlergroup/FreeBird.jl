# FreeBird.jl

Documentation for FreeBird.jl

## Installation

FreeBird.jl is a [Julia](http://julialang.org) package. Install `Julia` first following the instructions in the [download page](https://julialang.org/downloads/).

Once Julia is installed, you can install the `FreeBird` package from the Julia REPL with:

```julia-repl
julia> ] # press the "]" key on your keyboard to enter the Pkg manager
```
then

```julia-repl
pkg> add FreeBird
```
Alternatively, you can install the package with:
```julia
using Pkg; Pkg.add("FreeBird")
```
Useful for installing the package in a script or a notebook.

Or, if you want to install a specific branch from GitHub, you can do so with:

```julia-repl
pkg> add https://github.com/wexlergroup/FreeBird.jl#branch_name
```
Or, again, to install the package in a script or a notebook:
```julia
using Pkg; Pkg.add(url="https://github.com/wexlergroup/FreeBird.jl",rev="branch_name")
```
Remember to replace `branch_name` with the name of the branch you want to install.

To get back to the Julia REPL, press `Ctrl+C` or backspace (when the REPL cursor is at the beginning of the input).

