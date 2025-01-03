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

pkg> add FreeBird
```
Or, if you want to install a specific branch from GitHub, you can do so with:

```julia-repl
pkg> add https://github.com/wexlergroup/FreeBird.jl#branch_name
```
To get back to the Julia REPL, press `Ctrl+C` or backspace (when the REPL cursor is at the beginning of the input).

## References

Published surface nested sampling work from Wexler group:
- Yang, M.; B. Pártay, L.; B. Wexler, R. Surface Phase Diagrams from Nested Sampling. *Phys. Chem. Chem. Phys.* **2024**. https://doi.org/10.1039/D4CP00050A.

The nested sampling method:
- Pártay, L. B.; Csányi, G.; Bernstein, N. Nested Sampling for Materials. *Eur. Phys. J. B* **2021**, 94 (8), 159. https://doi.org/10.1140/epjb/s10051-021-00172-1.

- Ashton, G.; Bernstein, N.; Buchner, J.; Chen, X.; Csányi, G.; Fowlie, A.; Feroz, F.; Griffiths, M.; Handley, W.; Habeck, M.; Higson, E.; Hobson, M.; Lasenby, A.; Parkinson, D.; Pártay, L. B.; Pitkin, M.; Schneider, D.; Speagle, J. S.; South, L.; Veitch, J.; Wacker, P.; Wales, D. J.; Yallup, D. Nested Sampling for Physical Scientists. *Nat Rev Methods Primers* **2022**, 2 (1), 39. https://doi.org/10.1038/s43586-022-00121-x.