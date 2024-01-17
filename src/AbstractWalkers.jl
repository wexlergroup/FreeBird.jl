module AbstractWalkers

using ExtXYZ
using LinearAlgebra
using ..Potentials

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart
export free_free_energy, frozen_energy, free_frozen_energy
export interaction_energy, total_energy

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

function free_free_energy(at::Atoms, lj::LJParameters; frozen::Int64=0)
    free_free_energy = 0.0
    for i in (frozen+1):length(at)
        for j in (i+1):length(at)
            r = norm(at.atom_data.position[i] - at.atom_data.position[j]).val
            free_free_energy += lj_energy(r,lj)
        end
    end
    return free_free_energy
end

function frozen_energy(at::Atoms, lj::LJParameters, frozen::Int64)
    frozen_energy = 0.0
    for i in 1:frozen
        for j in (i+1):frozen
            r = norm(at.atom_data.position[i] - at.atom_data.position[j]).val
            frozen_energy += lj_energy(r,lj)
        end
    end
    return frozen_energy
end

function free_frozen_energy(at::Atoms, lj::LJParameters, frozen::Int64)
    free_frozen_energy = 0.0
    for i in 1:frozen
        for j in (frozen+1):length(at)
            r = norm(at.atom_data.position[i] - at.atom_data.position[j]).val
            free_frozen_energy += lj_energy(r,lj)
            # @debug println("free_frozen_energy: ", free_frozen_energy, " r: ", r, " i and j: ", i, " ", j)
        end
    end
    return free_frozen_energy
end

function interaction_energy(at::Atoms, lj::LJParameters; frozen::Int64=0)
    if frozen == 0
        return free_free_energy(at, lj)
    else
        return free_free_energy(at, lj; frozen=frozen) + free_frozen_energy(at, lj, frozen)
    end
end

function total_energy(at::Atoms, lj::LJParameters; frozen::Int64=0)
    if frozen == 0
        return free_free_energy(at, lj)
    else
        return interaction_energy(at, lj; frozen=frozen) + frozen_energy(at, lj, frozen)
    end
end



end # module AbstractWalkers