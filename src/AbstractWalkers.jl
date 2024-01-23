module AbstractWalkers

using Setfield
using ..Potentials
using ..EnergyEval

export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart

abstract type AtomWalkers end

function assign_lj_energies!(ats::Vector{T}, lj::LJParameters; frozen::Int64=0, e_frozen::Float64=0.0) where T
    for i in eachindex(ats)
        at = ats[i]
        e_total = interaction_energy(at, lj; frozen=frozen) + e_frozen
        k = keys(at.system_data)
        v = values(at.system_data)
        # check if energy is already in system_data
        if !hasproperty(at.system_data, :energy)
            @debug "energy not in system_data, adding it"
            k = (k..., :energy)
            v = (v..., e_total)
            at = @set at.system_data = NamedTuple{k}(v)
            ats = @set ats[i] = at
        else
            at = @set at.system_data.energy = e_total
            ats[i] = at
        end
        # check if iter is already in system_data
        if !hasproperty(at.system_data, :iter)
            @debug "iter not in system_data, adding it"
            k = (k..., :iter)
            v = (v..., 0::Int)
            at = @set at.system_data = NamedTuple{k}(v)
            ats = @set ats[i] = at
        else
            at = @set at.system_data.iter = 0
            ats[i] = at
        end
    end
    return ats, typeof(ats[1])
end

struct LJAtomWalkers{T} <: AtomWalkers
    walkers::Vector{T}
    lj_potential::LJParameters
    function LJAtomWalkers(walkers::Vector{T}, lj_potential::LJParameters) where T
        walkers, new_T = assign_lj_energies!(walkers, lj_potential)
        return new{new_T}(walkers, lj_potential)
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
    walkers, new_T = assign_lj_energies!(walkers, lj_potential; frozen=num_frozen_particles, e_frozen=energy_frozen_particles)
    return LJAtomWalkersWithFrozenPart{new_T}(walkers, lj_potential, num_frozen_particles, energy_frozen_particles)
end



end # module AbstractWalkers