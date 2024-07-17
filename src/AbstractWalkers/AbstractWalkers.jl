"""
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using Combinatorics
using LinearAlgebra
using Statistics
using ..Potentials
# using ..EnergyEval
using ..Hamiltonians

export AbstractWalker
export AtomWalker, AtomWalkers, LJAtomWalkers
export LatticeWalkers, LatticeSystem, LatticeWalker  # , LatticeWalkers
export update_walker!
export exact_enumeration, nvt_monte_carlo, wang_landau

abstract type AbstractWalker end

include("AtomWalkers.jl")

include("LatticeWalkers.jl")

end # module AbstractWalkers