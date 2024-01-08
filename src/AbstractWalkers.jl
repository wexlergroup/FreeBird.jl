module AbstractWalkers

# using ExtXYZ
using ..Potentials

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart

abstract type AtomWalkers end

struct LJAtomWalkers{T} <: AtomWalkers
    walkers::Vector{T}
    lj_potential::LJParameters
end

struct LJAtomWalkersWithFrozenPart{T} <: AtomWalkers
    walkers::Vector{T}
    lj_potential::LJParameters
    num_frozen_particles::Int64
    energy_frozen_particles::Float64
end

end # module AbstractWalkers