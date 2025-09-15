# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                       frozen_energy definitions                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# multi-component frozen energy:
"""
    frozen_energy(at::AbstractSystem, cps::CompositeParameterSets{C,P}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy of the frozen particles in the system using a multi-component potential (CompositeParameterSets).
I.e., the components interact with each other using different parameters.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `cps::CompositeParameterSets{C,P}`: The composite potential parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       cps::CompositeParameterSets{C,P}, 
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       ) where {C,P}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], cps.param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in 1:C
            if frozen[i] && frozen[j] && i != j # both frozen and different
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], cps.param_sets[i,j])
            end
        end
    end
    return energy
end

# single-component frozen energy:
"""
    frozen_energy(at::AbstractSystem, pot::SingleComponentPotential{S}, list_num_par::Vector{Int}, frozen::Vector{Bool}) where {S}

Calculate the energy of the frozen particles in the system using a single-component potential.
I.e., the components interact with each other using the same parameters.
The energy is calculated by summing the pairwise interactions between the frozen particles.
Since the frozen particles do not move, the energy is typically only calculated once for a given system.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `pot::SingleComponentPotential{S}`: The single-component potential.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy of the frozen particles in the system.

"""
function frozen_energy(at::AbstractSystem, 
                       pot::SingleComponentPotential{S},
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool}
                       ) where {S}
    if length(list_num_par) != length(frozen)
        throw(ArgumentError("The number of frozen and free parts does not match the length of the number of components."))
    end
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(frozen) # find frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], pot)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if frozen[i] && frozen[j] # both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], pot)
            end
        end
    end
    return energy
end



# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                  interacting_energy definitions                            ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# single-component interacting energy:

"""
    interacting_energy(at::AbstractSystem, pot::SingleComponentPotential{S})

Calculate the energy from interactions between particles using the a single-component potential.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `pot::SingleComponentPotential{S}`: The single-component potential.

# Returns
- `energy`: The energy from interactions between particles.

"""
interacting_energy(at::AbstractSystem, pot::SingleComponentPotential{S}) where S = intra_component_energy(at, pot)
