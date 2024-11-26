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
export SLattice, LatticeWalker
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export MLattice
export update_walker!
export num_sites, occupied_site_count

abstract type AbstractWalker end

include("AtomWalkers.jl")

include("LatticeWalkers.jl")

include("helpers.jl")

end # module AbstractWalkers