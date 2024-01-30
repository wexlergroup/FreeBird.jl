module AbstractWalkers

using Setfield
using AtomsBase
using ..Potentials
using ..EnergyEval

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart

abstract type AtomWalkers end

function assign_lj_energies!(ats::Vector{T}, lj::LJParameters; frozen::Int64=0, e_frozen::Float64=0.0) where T
    for i in eachindex(ats)
        e_total = interaction_energy(ats[i], lj; frozen=frozen) + e_frozen
        ats[i] = FlexibleSystem(ats[i].particles, ats[i].bounding_box, ats[i].boundary_conditions; iter = 0, energy = e_total)
    end
    return ats
end

struct LJAtomWalkers{T} <: AtomWalkers
    walkers::Vector{T}
    lj_potential::LJParameters
    function LJAtomWalkers(walkers::Vector{T}, lj_potential::LJParameters) where T
        assign_lj_energies!(walkers, lj_potential)
        return new{T}(walkers, lj_potential)
    end
end

struct LJAtomWalkersWithFrozenPart{T} <: AtomWalkers
    walkers::Vector{T}
    lj_potential::LJParameters
    num_frozen_particles::Int64
    energy_frozen_particles::Float64
end

function LJAtomWalkersWithFrozenPart(walkers::Vector{T}, lj_potential::LJParameters, num_frozen_particles::Int64) where T
    energy_frozen_particles = frozen_energy(walkers[1], lj_potential, num_frozen_particles)
    assign_lj_energies!(walkers, lj_potential; frozen=num_frozen_particles, e_frozen=energy_frozen_particles)
    return LJAtomWalkersWithFrozenPart{T}(walkers, lj_potential, num_frozen_particles, energy_frozen_particles)
end



end # module AbstractWalkers