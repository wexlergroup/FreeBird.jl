"""
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using LinearAlgebra
using Statistics
using ..Potentials
using ..Hamiltonians

export AbstractWalker
export AbstractLattice
export AtomWalker
export LatticeWalker
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export MLattice, SLattice, GLattice
export update_walker!
export num_sites, occupied_site_count

abstract type AbstractWalker end

include("atomistic_walkers.jl")

include("lattice_walkers.jl")

include("helpers.jl")

end # module AbstractWalkers