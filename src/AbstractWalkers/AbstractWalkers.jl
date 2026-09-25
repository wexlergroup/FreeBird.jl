"""
    AbstractWalkers
    
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using LinearAlgebra
using Statistics
using ASEconvert
using PythonCall
using ..AbstractPotentials
using ..AbstractHamiltonians

const _PY_COPY = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(_PY_COPY, pyimport("copy"))
end

export AbstractWalker
export AbstractLattice
export AtomWalker
export sort_components_by_atomic_number
export split_components
export split_components_by_chemical_species
export check_num_components
export insert_particle!, remove_particle!
export LatticeWalker
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export MLattice, SLattice, GLattice, AtomicLattice
export replicate_walkers
export update_walker!
export num_sites, occupied_site_count, num_lattice_components
export coverage, sync_ase_lattice!, nn_distance
export n_occupied, is_occupied, set_occupied!, occupied_indices, empty_indices
export swap_sites!, neighbor_shell
export order_parameter_c2x2
export order_parameter_sqrt3
export bragg_amplitude
export order_parameter_stripe
export order_parameter_p2x2
export bragg_amplitude_layers, stacking_bragg_amplitude
export site_layers, layer_coverage, occupancy_profile, layer_field
export enumerate_motif_embeddings, motif_distances
export view_structure

abstract type AbstractWalker end

include("atomistic_walkers.jl")

include("lattice_walkers.jl")

include("helpers.jl")

include("shows.jl")

include("lattice_interface.jl")

end # module AbstractWalkers
