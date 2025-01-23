"""
    AbstractLiveSets

Module for defining the livesets, which are collections of 
walkers that are used in the sampling schemes.
"""
module AbstractLiveSets

import Unitful: unit

using ..AbstractWalkers
using ..EnergyEval
using ..AbstractPotentials
using ..AbstractHamiltonians

export AbstractLiveSet
export AtomWalkers, LatticeWalkers
export LJAtomWalkers, LatticeGasWalkers

abstract type AbstractLiveSet end


include("atomistic_livesets.jl")

include("lattice_livesets.jl")
    
end # module AbstractLiveSets