"""
    MonteCarloMoves

Module containing functions for performing Monte Carlo moves on atomistic/lattice systems.
"""
module MonteCarloMoves

using ExtXYZ
using AtomsBase
using Distributions
using Unitful
using StaticArrays

using ..AbstractPotentials
using ..AbstractHamiltonians
using ..AbstractWalkers
using ..EnergyEval

export periodic_boundary_wrap!
export single_atom_random_walk!, two_atoms_swap!
export MC_random_walk!, MC_random_walk_2D!
export MC_new_sample!, MC_rejection_sampling!, MC_random_swap!
export lattice_random_walk!
export generate_random_new_lattice_sample!
export MC_mixed_moves!
export MC_cluster_walk!
export geometric_cluster_swap!
export random_microstate!
export lattice_insert_particle!, lattice_delete_particle!
export lattice_biased_sites
export MC_grand_canonical_walk!
export MC_muVT_walk!

export free_component_index, free_par_index

include("helpers.jl")

include("cluster_moves.jl")

include("random_walks.jl")

include("atomistic_swaps.jl")

include("mixed_moves.jl")

end # module MonteCarloMoves