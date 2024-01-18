module AbstractWalkers

using ExtXYZ
using LinearAlgebra
using Setfield
using ..Potentials
using ..EnergyEval

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart
export assign_lj_energies!

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

function LJAtomWalkersWithFrozenPart(walkers::Vector{T}, lj_potential::LJParameters, num_frozen_particles::Int64) where T
    energy_frozen_particles = frozen_energy(walkers[1], lj_potential, num_frozen_particles)
    return LJAtomWalkersWithFrozenPart{T}(walkers, lj_potential, num_frozen_particles, energy_frozen_particles)
end



"""
    assign_lj_energies!(liveset::LJAtomWalkers)

Assigns the Lennard-Jones potential energies to the walkers in the given `liveset`.

# Arguments
- `liveset::LJAtomWalkers`: The liveset of walkers to be assigned energies.
"""
function assign_lj_energies!(liveset::LJAtomWalkers)
    ats = liveset.walkers
    lj = liveset.lj_potential
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = total_energy(at, lj)
        ats[i] = at
    end
end

function assign_lj_energies!(liveset::LJAtomWalkersWithFrozenPart)
    ats = liveset.walkers
    lj = liveset.lj_potential
    frozen = liveset.num_frozen_particles
    e_frozen = liveset.energy_frozen_particles
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = interaction_energy(at, lj; frozen=frozen) + e_frozen
        ats[i] = at
    end
end



end # module AbstractWalkers