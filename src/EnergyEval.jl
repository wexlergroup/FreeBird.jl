"""
Module for evaluating energy-related quantities for a system.
"""
module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using ..Potentials

export free_free_energy, frozen_energy, free_frozen_energy
export interaction_energy, total_energy
export pbc_dist

"""
    pbc_dist(pos1, pos2, at)

Compute the distance between two positions considering periodic boundary conditions.

# Arguments
- `pos1::Union{SVector{T},Vector{T}}`: The first position.
- `pos2::Union{SVector{T},Vector{T}}`: The second position.
- `at::AbstractSystem`: The abstract system containing boundary conditions and bounding box.

# Returns
- `d::Float64`: The distance between `pos1` and `pos2` considering periodic boundary conditions.

"""
function pbc_dist(pos1::Union{SVector{T},Vector{T}},
                      pos2::Union{SVector{T},Vector{T}},
                      at::AbstractSystem) where {T}
    pbc = at.boundary_conditions
    box = at.bounding_box
    distsq = 0.0u"Å"^2
    for i in eachindex(pos1)
        if pbc[i] == Periodic()
            distsq += min(abs(pos1[i] - pos2[i]), box[i][i] - abs(pos1[i] - pos2[i]))^2
        elseif pbc[i] == DirichletZero()
            distsq += (pos1[i] - pos2[i])^2
        else
            error("Unsupported boundary condition: $(pbc[i])")
        end
    end
    return sqrt(distsq)
end

"""
    free_free_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)

Calculate the energy from interactions between free particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `frozen::Int64`: The number of frozen particles in the system.

# Returns
- `free_free_energy`: The energy from interactions between free particles.

"""
function free_free_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)
    free_free_energy = 0.0u"eV"
    for i in (frozen+1):length(at)
        for j in (i+1):length(at)
            r = pbc_dist(position(at, i), position(at, j), at)
            free_free_energy += lj_energy(r,lj)
        end
    end
    return free_free_energy
end

"""
    frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)

Compute the energy of the frozen particles in the system using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system containing the particles.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `frozen::Int64`: The number of frozen particles.

# Returns
- `frozen_energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)
    frozen_energy = 0.0u"eV"
    for i in 1:frozen
        for j in (i+1):frozen
            r = pbc_dist(position(at, i), position(at, j), at)
            frozen_energy += lj_energy(r,lj)
        end
    end
    return frozen_energy
end

"""
    free_frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)

Compute the energy from interactions between free and frozen particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free and frozen particles.

# Arguments
- `at::AbstractSystem`: The system containing the particles.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `frozen::Int64`: The number of frozen particles.

# Returns
- `free_frozen_energy`: The energy from interactions between free and frozen particles.

"""
function free_frozen_energy(at::AbstractSystem, lj::LJParameters, frozen::Int64)
    free_frozen_energy = 0.0u"eV"
    for i in 1:frozen
        for j in (frozen+1):length(at)
            r = pbc_dist(position(at, i), position(at, j), at)
            free_frozen_energy += lj_energy(r,lj)
        end
    end
    return free_frozen_energy
end

"""
    interaction_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)

Compute the energy from interactions between particles using the Lennard-Jones potential.
The energy is calculated by combining the energies from interactions between free particles 
and between free and frozen particles. The energy from interactions between frozen particles
is not included, as it is typically treated as a constant energy offset. If the number of frozen
particles is zero, the energy calculated is equivalent to the total energy of the system.

# Arguments
- `at::AbstractSystem`: The system containing the particles.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `frozen::Int64`: The number of frozen particles.

# Returns
- `interaction_energy`: The energy from interactions between particles.

"""
function interaction_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)
    e_free_frozen = frozen == 0 ? 0.0u"eV" : free_frozen_energy(at, lj, frozen)
    return free_free_energy(at, lj; frozen=frozen) + e_free_frozen
end


    
end # module EnergyEval