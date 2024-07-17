module AbstractLiveSets

using ..AbstractWalkers
using ..EnergyEval
using ..Potentials
using ..Hamiltonians

export AbstractLiveSet
export AtomWalkers, LatticeWalkers

abstract type AbstractLiveSet end

# Base.eltype(::AbstractLiveSet{T}) where {T} = T

include("AtomLiveSets.jl")

include("LatticeLiveSets.jl")
    
end # module AbstractLiveSets