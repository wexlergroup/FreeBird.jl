"""
    AbstractLiveSets

Module for defining the livesets, which are collections of 
walkers that are used in the sampling schemes.
"""
module AbstractLiveSets

using Distributed

import Unitful: unit

using ..AbstractWalkers
using ..EnergyEval
using ..AbstractPotentials
using ..AbstractHamiltonians

export AbstractLiveSet
export AtomWalkers, LatticeWalkers
export LJAtomWalkers, LatticeGasWalkers
export LJSurfaceWalkers
export MLIPAtomWalkers
export GuptaAtomWalkers

abstract type AbstractLiveSet end


include("atomistic_livesets.jl")

include("lattice_livesets.jl")

include("shows.jl")
    
end # module AbstractLiveSets