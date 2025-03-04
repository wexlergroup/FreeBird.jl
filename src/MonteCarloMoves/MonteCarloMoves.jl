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
export single_atom_random_walk!
export MC_random_walk!, MC_random_walk_2D!
export MC_new_sample!, MC_rejection_sampling!, MC_random_swap!
export generate_random_new_lattice_sample!

include("helpers.jl")

include("random_walks.jl")

include("atomistic_swaps.jl")

end # module MonteCarloMoves