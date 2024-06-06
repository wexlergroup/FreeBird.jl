"""
Module for evaluating energy-related quantities for a system.
"""
module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using ..Potentials

export pbc_dist
export interacting_energy, frozen_energy

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
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2) if i <= j]
    for (i, j) in pairs
        # @show i,j # DEBUG
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy += lj_energy(r,lj)
    end
    return energy
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
    split_components(at::AbstractSystem, list_num_par::Vector{Int})

Split the system into components based on the number of particles in each component.

# Arguments
- `at::AbstractSystem`: The system to split.
- `list_num_par::Vector{Int}`: The number of particles in each component.

# Returns
- `components`: An array of `FastSystem` objects representing the components of the system.

"""
function split_components(at::AbstractSystem, list_num_par::Vector{Int})
    components = Array{FastSystem}(undef, length(list_num_par))
    comp_cut = vcat([0],cumsum(list_num_par))
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    for i in 1:length(list_num_par)
        components[i] = FastSystem(at[comp_split[i]],at.bounding_box,at.boundary_conditions)
    end
    return components
end

"""
    check_num_components(C::Int, list_num_par::Vector{Int}, frozen::Vector{Bool})

Check that the number of components matches the length of the list of number of particles and frozen particles.

# Arguments
- `C::Int`: The number of components.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

"""
function check_num_components(C::Int, list_num_par::Vector{Int}, frozen::Vector{Bool})
    if length(list_num_par) != C
        throw(ArgumentError("The number of components does not match the length of the list of number of particles."))
    elseif length(frozen) != C
        throw(ArgumentError("The number of components does not match the length of the list of frozen particles."))
    end
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
        for j in (i+1):C
            if frozen[i] && frozen[j] # both frozen
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

end # module EnergyEval