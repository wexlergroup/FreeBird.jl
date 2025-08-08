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
using CUDA
export pbc_dist
export interacting_energy, frozen_energy
export single_site_energy
export extract_gpu_params
export LJParams
export pair_energy_gpu
export pbc_dist_gpu
export intra_component_energy
export compute_intra_pair_energies_kernel!
export LJParams


include("atomistic_energies.jl")

include("lattice_energies.jl")

end # module EnergyEval