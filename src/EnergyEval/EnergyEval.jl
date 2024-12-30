"""
Module for evaluating energy-related quantities for a system.
"""
module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using ..AbstractWalkers
using ..AbstractPotentials
using ..AbstractHamiltonians

export pbc_dist
export interacting_energy, frozen_energy
export sort_components_by_atomic_number
export split_components

include("atomistic_energies.jl")

include("lattice_energies.jl")

include("helpers.jl")

end # module EnergyEval