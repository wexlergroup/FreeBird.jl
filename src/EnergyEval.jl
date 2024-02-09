module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using ..Potentials

export free_free_energy, frozen_energy, free_frozen_energy
export interaction_energy, total_energy
export pbc_dist

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

function interaction_energy(at::AbstractSystem, lj::LJParameters; frozen::Int64=0)
    e_free_frozen = frozen == 0 ? 0.0u"eV" : free_frozen_energy(at, lj, frozen)
    return free_free_energy(at, lj; frozen=frozen) + e_free_frozen
end


    
end # module EnergyEval