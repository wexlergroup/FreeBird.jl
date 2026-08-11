"""
    struct IdealGasParameters

A zero-interaction potential representing a classical ideal gas.

`IdealGasParameters` is a `SingleComponentPotential{Pairwise}` whose `pair_energy`
is identically zero at any separation. All existing pairwise dispatch paths
(`intra_component_energy`, `inter_component_energy`, `interacting_energy`,
`frozen_energy`) therefore evaluate to `0.0u"eV"` without any further methods.

The intended use is as an analytical reference for grand-canonical post-processing:
the ideal-gas grand partition function has the closed form `Ξ(μ, T) = exp(zV)` with
`z = exp(βμ)/Λ(T)³`, which provides the cleanest possible check of normalization,
fugacity, and thermal-wavelength conventions in downstream analysis code.
"""
struct IdealGasParameters <: SingleComponentPotential{Pairwise} end

"""
    pair_energy(r::typeof(1.0u"Å"), ::IdealGasParameters)

Return `0.0u"eV"` at any separation `r`. The ideal gas has no pair interactions.
"""
pair_energy(::typeof(1.0u"Å"), ::IdealGasParameters) = 0.0u"eV"
