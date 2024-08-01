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
export AtomWalker, AtomWalkers, LJAtomWalkers
export LatticeWalkers, LatticeSystem, LatticeWalker  # , LatticeWalkers
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export update_walker!

abstract type AbstractWalker end

include("AtomWalkers.jl")

include("LatticeWalkers.jl")

include("helpers.jl")

end # module AbstractWalkers