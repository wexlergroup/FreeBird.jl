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
using ..AbstractPotentials
using ..AbstractHamiltonians

export AbstractWalker
export AbstractLattice
export AtomWalker
export sort_components_by_atomic_number
export split_components
export split_components_by_chemical_species
export check_num_components
export LatticeWalker
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export MLattice, SLattice, GLattice
export update_walker!
export num_sites, occupied_site_count
export view_structure

abstract type AbstractWalker end

include("atomistic_walkers.jl")

include("lattice_walkers.jl")

include("helpers.jl")

include("shows.jl")

end # module AbstractWalkers