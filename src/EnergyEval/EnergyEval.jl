"""
Module for evaluating energy-related quantities for a system.
"""
module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using ..AbstractWalkers
using ..Potentials
using ..Hamiltonians

export pbc_dist
export interacting_energy, frozen_energy
export sort_components_by_atomic_number
export split_components

include("helpers.jl")

"""
    pbc_dist(pos1, pos2, at)

Compute the distance between two positions considering periodic boundary conditions. Currently only works for orthorhombic lattices.

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
    inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, lj::LJParameters)

Compute the energy between two components of a system using the Lennard-Jones potential.

# Arguments
- `at1::AbstractSystem`: The first component of the system.
- `at2::AbstractSystem`: The second component of the system.
- `lj::LJParameters`: The Lennard-Jones parameters.

# Returns
- `energy`: The energy between the two components.

"""
function inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, lj::LJParameters)
    energy = 0.0u"eV"
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    for (i, j) in pairs
        # @show i,j # DEBUG
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy += lj_energy(r,lj)
    end
    return energy
end

function inter_component_energy_threaded(at1::AbstractSystem, at2::AbstractSystem, lj::LJParameters)
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    energy = Array{typeof(0.0u"eV"), 1}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        # @show i,j # DEBUG
        (i, j) = pairs[k]
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy[k] = lj_energy(r,lj)
    end
    # energy = energy*u"eV"
    return sum(energy)
end

"""
    intra_component_energy(at::AbstractSystem, lj::LJParameters)

Compute the energy within a component of a system using the Lennard-Jones potential.

# Arguments
- `at::AbstractSystem`: The component of the system.
- `lj::LJParameters`: The Lennard-Jones parameters.

# Returns
- `energy`: The energy within the component.

"""
function intra_component_energy(at::AbstractSystem, lj::LJParameters)
    energy = 0.0u"eV"
    for i in 1:length(at)
        for j in (i+1):length(at)
            r = pbc_dist(position(at, i), position(at, j), at)
            energy += lj_energy(r,lj)
        end
    end
    return energy
end



"""
    frozen_energy(at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy of the frozen particles in the system using a composite Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       ljs::CompositeLJParameters{C}, 
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       ) where {C}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], ljs.lj_param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in 1:C
            if frozen[i] && frozen[j] && i != j # both frozen and different
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], ljs.lj_param_sets[i,j])
            end
        end
    end
    return energy
end

"""
    frozen_energy(at::AbstractSystem, lj::LJParameters, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy of the frozen particles in the system using a single Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       lj::LJParameters,
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       )
    if length(list_num_par) != length(frozen)
        throw(ArgumentError("The number of frozen and free parts does not match the length of the number of components."))
    end
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], lj)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if frozen[i] && frozen[j] # both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], lj)
            end
        end
    end
    return energy
end


"""
    interacting_energy(at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            ljs::CompositeLJParameters{C}, 
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool}
                            ) where {C}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], ljs.lj_param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in (i+1):C
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], ljs.lj_param_sets[i,j])
            end
        end
    end
    return energy
end

"""
    interacting_energy(at::AbstractSystem, lj::LJParameters, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            lj::LJParameters,
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool}
                            )
    if length(list_num_par) != length(frozen)
        throw(ArgumentError("The number of frozen and free parts does not match the length of the number of components."))
    end
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], lj)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], lj)
            end
        end
    end
    return energy
end

"""
    interacting_energy(at::AbstractSystem, lj::LJParameters)

Calculate the energy from interactions between particles using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.

# Returns
- `energy`: The energy from interactions between particles.

"""
interacting_energy(at::AbstractSystem, lj::LJParameters) = intra_component_energy(at, lj)

"""
    interacting_energy(lattice::SLattice, h::LatticeGasHamiltonian)

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice::SLattice`: The lattice configuration.
- `h::LatticeGasHamiltonian`: The lattice-gas Hamiltonian parameters.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""
function interacting_energy(lattice::SLattice, h::LatticeGasHamiltonian)
    e_adsorption = sum(lattice.occupations .& lattice.adsorptions) * h.adsorption_energy
    e_nn = 0.0u"eV"
    e_nnn = 0.0u"eV"

    for index in 1:length(lattice.occupations)
        if lattice.occupations[index]
            # Compute nearest-neighbor interaction energy
            for nn in lattice.neighbors[index][1]
                if lattice.occupations[nn]
                    e_nn += h.nn_interaction_energy / 2
                end
            end

            # Compute next-nearest-neighbor interaction energy
            for nnn in lattice.neighbors[index][2]
                if lattice.occupations[nnn]
                    e_nnn += h.nnn_interaction_energy / 2
                end
            end
        end
    end

    e_interaction = e_adsorption + e_nn + e_nnn
    return e_interaction
end

"""
    lattice_interaction_energy(lattice_occupations::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U})

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice_occupations::Vector{Bool}`: The lattice occupation configuration.
- `lattice_neighbors::Vector{Vector{Vector{Int64}}}`: The lattice neighbor list.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::U`: The interaction energy of the lattice configuration.

"""
function lattice_interaction_energy(lattice_occupations::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    e_interaction::U = 0.0*unit(h.on_site_interaction)
    for index in eachindex(lattice_occupations)
        if lattice_occupations[index]
            for n in 1:N
                for neighbor in lattice_neighbors[index][n]
                    if lattice_occupations[neighbor]
                        e_interaction += h.nth_neighbor_interactions[n] / 2
                    end
                end
            end
        end
    end
    return e_interaction
end

"""
    inter_component_energy(lattice1::Vector{Bool}, lattice2::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U})

Compute the interaction energy between two lattice configurations using the Hamiltonian parameters.

# Arguments
- `lattice1::Vector{Bool}`: The first lattice configuration.
- `lattice2::Vector{Bool}`: The second lattice configuration.
- `lattice_neighbors::Vector{Vector{Vector{Int64}}}`: The lattice neighbor list.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::U`: The interaction energy between the two lattice configurations.

"""
function inter_component_energy(lattice1::Vector{Bool}, lattice2::Vector{Bool}, lattice_neighbors::Vector{Vector{Vector{Int64}}}, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    e_interaction::U = 0.0*unit(h.on_site_interaction)
    for index in eachindex(lattice1)
        if lattice1[index]
            for n in 1:N
                for neighbor in lattice_neighbors[index][n]
                    if lattice2[neighbor]
                        e_interaction += h.nth_neighbor_interactions[n]
                    end
                end
            end
        end
    end
    return e_interaction
end

"""
    interacting_energy(lattice::SLattice, h::GenericLatticeHamiltonian{N})

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice::SLattice`: The lattice configuration.
- `h::GenericLatticeHamiltonian{N,U}`: The generic lattice Hamiltonian parameters.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""
function interacting_energy(lattice::SLattice, h::GenericLatticeHamiltonian{N,U}) where {N,U}
    e_interaction::U = lattice_interaction_energy(lattice.components[1], lattice.neighbors, h)
    e_adsorption::U = sum(lattice.components[1] .& lattice.adsorptions) * h.on_site_interaction
    return e_interaction + e_adsorption
end

function interacting_energy(lattice::SLattice, h::MLatticeHamiltonian{C,N,U}) where {C,N,U}
    ham = h.Hamiltonians[1,1]
    e_interaction::U = lattice_interaction_energy(lattice.components[1], lattice.neighbors, ham)
    e_adsorption::U = sum(lattice.components[1] .& lattice.adsorptions) * ham.on_site_interaction
    return e_interaction + e_adsorption
end

"""
    interacting_energy(lattice::MLattice{C,G}, h::MLatticeHamiltonian{C,N,U})

Compute the interaction energy of a multi-component lattice configuration using the Hamiltonian parameters.

# Arguments
- `lattice::MLattice{C,G}`: The multi-component lattice configuration.
- `h::MLatticeHamiltonian{C,N,U}`: The multi-component lattice Hamiltonian parameters.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""
function interacting_energy(lattice::MLattice{C,G}, h::MLatticeHamiltonian{C,N,U}) where {C,G,N,U}
    adsorption_energy = 0.0*unit(h.Hamiltonians[1,1].on_site_interaction)
    interaction_energy = 0.0*unit(h.Hamiltonians[1,1].on_site_interaction)
    for i in 1:C
        interaction_energy += lattice_interaction_energy(lattice.components[i], lattice.neighbors, h.Hamiltonians[i,i])
        adsorption_energy += sum(lattice.components[i] .& lattice.adsorptions) * h.Hamiltonians[i,i].on_site_interaction
        for j in (i+1):C
            interaction_energy += inter_component_energy(lattice.components[i], lattice.components[j], lattice.neighbors, h.Hamiltonians[i,j])
        end
    end
    @debug "interaction energy: $interaction_energy, adsorption energy: $adsorption_energy"
    return interaction_energy + adsorption_energy
end

end # module EnergyEval