module AbstractLiveSets

import Unitful: unit

using ..AbstractWalkers
using ..EnergyEval
using ..Potentials
using ..Hamiltonians

export AbstractLiveSet
export AtomWalkers, LatticeWalkers
export LJAtomWalkers, LatticeGasWalkers

abstract type AbstractLiveSet end

# Base.eltype(::AbstractLiveSet{T}) where {T} = T

include("AtomLiveSets.jl")

include("LatticeLiveSets.jl")
    
end # module AbstractLiveSets