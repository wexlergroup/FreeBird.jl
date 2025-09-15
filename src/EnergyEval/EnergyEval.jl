"""
    EnergyEval
    
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
export single_site_energy

# definitions of frozen_energy and interacting_energy
include("atomistic_energies.jl")

# definitions of pbc_dist for calculating distances with periodic boundary conditions
include("atomistic_pbc_dist.jl")

# definitions of interacting_energy for pairwise potentials
include("atomistic_pairwise.jl")

# definitions of single_site_energy for computing a site energy using a pairwise potential
include("atomistic_single_site.jl")

# definitions of interacting_energy for many-body potentials
include("atomistic_many_body.jl")

# definitions of interacting_energy for lattice systems
include("lattice_energies.jl")

end # module EnergyEval