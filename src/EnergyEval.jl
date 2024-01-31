module EnergyEval

using AtomsBase
using LinearAlgebra
using ..Potentials

export free_free_energy, frozen_energy, free_frozen_energy
export interaction_energy, total_energy

function free_free_energy(at::FlexibleSystem, lj::LJParameters; frozen::Int64=0)
    free_free_energy = 0.0
    for i in (frozen+1):length(at)
        for j in (i+1):length(at)
            r = norm(at.particles[i][:position] - at.particles[j][:position]).val
            free_free_energy += lj_energy(r,lj)
        end
    end
    return free_free_energy
end

function frozen_energy(at::FlexibleSystem, lj::LJParameters, frozen::Int64)
    frozen_energy = 0.0
    for i in 1:frozen
        for j in (i+1):frozen
            r = norm(at.particles[i][:position] - at.particles[j][:position]).val
            frozen_energy += lj_energy(r,lj)
        end
    end
    return frozen_energy
end

function free_frozen_energy(at::FlexibleSystem, lj::LJParameters, frozen::Int64)
    free_frozen_energy = 0.0
    for i in 1:frozen
        for j in (frozen+1):length(at)
            r = norm(at.particles[i][:position] - at.particles[j][:position]).val
            free_frozen_energy += lj_energy(r,lj)
            # @debug println("free_frozen_energy: ", free_frozen_energy, " r: ", r, " i and j: ", i, " ", j)
        end
    end
    return free_frozen_energy
end

function interaction_energy(at::FlexibleSystem, lj::LJParameters; frozen::Int64=0)
    if frozen == 0
        return free_free_energy(at, lj)
    else
        return free_free_energy(at, lj; frozen=frozen) + free_frozen_energy(at, lj, frozen)
    end
end

function total_energy(at::FlexibleSystem, lj::LJParameters; frozen::Int64=0)
    if frozen == 0
        return free_free_energy(at, lj)
    else
        return interaction_energy(at, lj; frozen=frozen) + frozen_energy(at, lj, frozen)
    end
end


    
end # module EnergyEval