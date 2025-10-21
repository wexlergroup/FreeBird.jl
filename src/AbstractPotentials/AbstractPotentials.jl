"""
    AbstractPotentials

Module for defining and implementing potentials.

Currently implemented potentials:
- Lennard-Jones potential
- Gupta (SMA-TB) potential
"""
module AbstractPotentials

using Unitful

export AbstractPotential
export SingleComponentPotential, MultiComponentPotential

export PotentialStyle
export Pairwise, ManyBody


export CompositeParameterSets

export LJParameters, lj_energy
export LennardJonesParameterSets
export GuptaParameters
export PyCalculator


using AtomsBase, StaticArrays, Unitful
using AtomsCalculators
using ASEconvert
using PythonCall

export ASELennardJones, MACEPotential, OrbPotential, MLPotential

export pair_energy, two_body_energy, many_body_energy, total_energy

""" 
    AbstractPotential
An abstract type for potentials. All potentials should be subtypes of `AbstractPotential`.

Currently, there are two main categories of potentials:  

    1. `SingleComponentPotential{T}`: Potentials that use the same parameters for all components. 
       Here, `T` is the type of the single component parameter sets, e.g., `LJParameters` or `GuptaParameters`.  

    2. `MultiComponentPotential`: Potentials that use different parameters for different components, e.g., `CompositeParameterSets{C,P}`.
       This type is not parameterized by `C` or `P` to allow for more flexibility in defining multi-component potentials, i.e.,
       using different interaction models for different pairs of components.
"""
abstract type AbstractPotential end

"""
    SingleComponentPotential{T}
An abstract type for single-component potentials, where all components interact using the same parameters.
Here, `T` is the type of the single component parameter sets, e.g., `LJParameters` or `GuptaParameters`.
"""
abstract type SingleComponentPotential{T} <: AbstractPotential end

"""
    MultiComponentPotential
An abstract type for multi-component potentials, where different components can interact using different parameters.
"""
abstract type MultiComponentPotential <: AbstractPotential end

"""
    PotentialStyle
An abstract type for potential styles, used to parameterize potentials and dispatch methods for energy evaluation based on the interaction model.

# Subtypes
- `Pairwise`: Potentials that depend only on pairwise interactions between particles, e.g., Lennard-Jones potential.
- `ManyBody`: Potentials that depend on many-body interactions, e.g., Gupta potential.

# Examples
- `Pairwise`: `LJParameters <: LennardJonesParameterSets{Pairwise}`.
- `ManyBody`: `GuptaParameters <: GuptaParameters{ManyBody}`.
"""
abstract type PotentialStyle end

"""
    Pairwise
A subtype of `PotentialStyle` representing pairwise interaction models.
"""
abstract type Pairwise <: PotentialStyle end

"""
    ManyBody
A subtype of `PotentialStyle` representing many-body interaction models.
"""
abstract type ManyBody <: PotentialStyle end

# defining how to construct CompositeParameterSets
include("composite_paramsets.jl")

# collection of interatomic potentials
# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                          Interatomic Potentials                            ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ - Lennard-Jones potential                                                  ║
# ║ - Gupta potential                                                          ║
# ╚════════════════════════════════════════════════════════════════════════════╝
include("lennardjones.jl")

include("gupta.jl")

include("ase_calculators.jl")

end # module AbstractPotentials