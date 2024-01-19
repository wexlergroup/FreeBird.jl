module AbstractWalkers

using Setfield
using ..Potentials
using ..EnergyEval

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart

abstract type AtomWalkers end

function assign_lj_energies!(ats::Vector{T}, lj::LJParameters) where T
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = total_energy(at, lj)
        ats[i] = at
    end
end

function assign_lj_energies!(ats::Vector{T}, lj::LJParameters, frozen::Int64, e_frozen::Float64) where T
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = interaction_energy(at, lj; frozen=frozen) + e_frozen
        ats[i] = at
    end
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
    assign_lj_energies!(walkers, lj_potential, num_frozen_particles, energy_frozen_particles)
    return LJAtomWalkersWithFrozenPart{T}(walkers, lj_potential, num_frozen_particles, energy_frozen_particles)
end



end # module AbstractWalkers