# FreeBird

[![Build Status](https://github.com/wexlergroup/FreeBird.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wexlergroup/FreeBird.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/wexlergroup/FreeBird.jl/badge.svg?branch=main)](https://coveralls.io/github/wexlergroup/FreeBird.jl?branch=main)
[![Dev](https://img.shields.io/badge/docs-stable-gre.svg)](https://wexlergroup.github.io/FreeBird.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wexlergroup.github.io/FreeBird.jl/dev/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

*Free energy calculators by Bayesian-inspired nested sampling and other integration techniques*

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

## References

The FreeBird.jl code paper in the Journal of Chemical Theory and Computation:
- Yang, R.; Chen, J.; Thibodeaux, D.; Wexler, R. B. FreeBird.jl: An Extensible Toolbox for Simulating Interfacial Phase Equilibria. *J. Chem. Theory Comput.* **2025**, *21* (21), 10765–10779. [https://doi.org/10.1021/acs.jctc.5c01348](https://doi.org/10.1021/acs.jctc.5c01348)

Published surface nested sampling work from the Wexler group:
- Yang, M.; Pártay, L. B.; Wexler, R. B. Surface phase diagrams from nested sampling. *Phys. Chem. Chem. Phys.* **2024**, *26* (18), 13862–13874. [https://doi.org/10.1039/D4CP00050A](https://doi.org/10.1039/D4CP00050A)
- Chatbipho, T.; Yang, R.; Wexler, R. B.; Pártay, Livia B. Adsorbate phase transitions on nanoclusters from nested sampling. *J. Chem. Phys.* **2025**, *163* (17), 174701. [https://doi.org/10.1063/5.0283538](https://doi.org/10.1063/5.0283538)

The nested sampling method:
- Pártay, L. B.; Bartók, A. P.; Csányi, G. Efficient Sampling of Atomic Configurational Spaces. *J. Phys. Chem. B* **2010**, *114* (32), 10502–10512. [https://doi.org/10.1021/jp1012973](https://doi.org/10.1021/jp1012973)
- Ashton, G.; Bernstein, N.; Buchner, J.; Chen, X.; Csányi, G.; Fowlie, A.; Feroz, F.; Griffiths, M.; Handley, W.; Habeck, M.; Higson, E.; Hobson, M.; Lasenby, A.; Parkinson, D.; Pártay, L. B.; Pitkin, M.; Schneider, D.; Speagle, J. S.; South, L.; Veitch, J.; Wacker, P.; Wales, D. J.; Yallup, D. Nested Sampling for Physical Scientists. *Nat, Rev. Methods Primers* **2022**, *2* (1), 39. [https://doi.org/10.1038/s43586-022-00121-x](https://doi.org/10.1038/s43586-022-00121-x)
