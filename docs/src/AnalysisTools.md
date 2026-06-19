# Analysis Tools

The `AnalysisTools` module post-processes sampling output into thermodynamic
quantities: phase-space weights (`ωᵢ`), the partition function, internal energy,
constant-volume heat capacity (`cv`), grand-canonical reweighting, and — described
below — microcanonical inflection-point analysis of phase transitions.

## Microcanonical inflection-point analysis

Nested sampling yields the microcanonical entropy ``S(E) = \ln g(E)`` almost for
free. Along the contiguous cull index ``i`` the enclosed prior volume is
``X_i = (K/(K+1))^i``, so ``\ln X_i = i\,\ln(K/(K+1))`` and, with the recorded
energy ladder ``E(i)``,

```math
S(E) = i\,\ln\!\frac{K}{K+1} - \ln\left|\frac{dE}{di}\right| ,
```

where derivatives are taken against the dense, uniform ``i``-axis (local cubic
fits) rather than a binned ``\ln g`` in energy space. Following Schnabel *et al.*,
Phys. Rev. E **84**, 011127 (2011) and Qi & Bachmann, Phys. Rev. Lett. **120**,
180601 (2018), phase transitions in finite systems appear as *inflection points* of
the inverse caloric temperature ``\beta(E) = dS/dE``, and their **order** follows
from the higher derivatives:

| transition order | signature in ``S^{(n)}(E)`` |
|:--|:--|
| 1 (first-order)  | ``\beta = S'`` has a positive local **minimum** (``\beta`` backbends; latent heat) |
| 2 (second-order) | ``\gamma = S''`` has a negative local **maximum** |
| 3                | ``\delta = S'''`` has a positive local minimum |
| odd ``n``        | positive local minimum of ``S^{(n)}`` |
| even ``n``       | negative local maximum of ``S^{(n)}`` |

The transition temperature is ``T_\mathrm{tr} = 1/\beta(E_\mathrm{tr})``. This is a
purely microcanonical route: it needs no temperature grid, returns transition
**energies and orders** directly, and can resolve transitions that broad canonical
`cv(T)` peaks blur.

### Workflow

```julia
using FreeBird

df = read_output("ladder.csv")            # canonical NS output (columns :iter, :emax)
K  = 320                                   # number of live walkers used to produce it

E, S = microcanonical_entropy(df, K)               # S(E) = ln g(E)
d    = caloric_derivatives(df, K; max_order = 2)   # d.E, d.β, d.γ
ts   = inflection_transitions(df, K; max_order = 2) # → (E_tr, T_tr, order, kind, strength)

# Recommended: confirm each transition is converged in the walker count K
dfs = [read_output("ladder_K320.csv"),
       read_output("ladder_K640.csv"),
       read_output("ladder_K1280.csv")]
conv = transition_convergence(dfs, [320, 640, 1280])   # flags drifting transitions
```

### Convergence and smoothing

!!! warning "Use enough walkers, and check convergence"
    The derivative ``\gamma = d^2S/dE^2`` amplifies the finite-walker "staircase"
    granularity of the NS ladder (step ``\sim 1/(K\,g)``). This roughness is
    *deterministic* in the walker count ``K`` — it does **not** average out over
    independent runs — so a transition should only be trusted once its temperature
    and order stop drifting as ``K`` grows; use [`transition_convergence`](@ref).
    Near-ground transitions converge last. The light internal smoothing
    (`halfwidth`) is a necessary differentiation aid, **not** a substitute for
    adequate ``K``: smoothing heavily to mask too-few walkers biases the transition
    temperatures. Energy-window/`ground_trim` controls remove the ``\beta``
    divergence as ``E \to E_\mathrm{ground}``.

## Functions
```@autodocs
Modules = [AnalysisTools]
```
