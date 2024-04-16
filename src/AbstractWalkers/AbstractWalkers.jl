"""
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using Combinatorics
using ..Potentials
using ..EnergyEval

export AtomWalker, AtomWalkers, LJAtomWalkers
export LatticeWalkers, Lattice2DWalker, Lattice2DWalkers
export update_walker!

include("AtomWalkers.jl")

include("LatticeWalkers.jl")

end # module AbstractWalkers