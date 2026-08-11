"""
    inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::SingleComponentPotential{Pairwise})

Compute the energy between two components of a system using a specified pairwise potential.

# Arguments
- `at1::AbstractSystem`: The first component of the system.
- `at2::AbstractSystem`: The second component of the system.
- `pot::SingleComponentPotential{Pairwise}`: The potential used to compute the energy.

# Returns
- `energy`: The energy between the two components.

"""
function inter_component_energy(at1::AbstractSystem, at2::AbstractSystem, pot::SingleComponentPotential{Pairwise})
    # build pairs of particles
    pairs = [(i, j) for i in 1:length(at1), j in 1:length(at2)]
    # @show pairs # DEBUG
    energy = Array{typeof(0.0u"eV"), 1}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        # @show i,j # DEBUG
        (i, j) = pairs[k]
        r = pbc_dist(position(at1, i), position(at2, j), at1)
        energy[k] = pair_energy(r, pot)
    end
    # energy = energy*u"eV"
    return sum(energy)
end

"""
    intra_component_energy(at::AbstractSystem, pot::SingleComponentPotential{Pairwise})

Compute the energy within a component of a system using a specified pairwise potential.

# Arguments
- `at::AbstractSystem`: The component of the system.
- `pot::SingleComponentPotential{Pairwise}`: The potential used to compute the energy.

# Returns
- `energy`: The energy within the component.

"""
function intra_component_energy(at::AbstractSystem, pot::SingleComponentPotential{Pairwise}) 
    # num_pairs = length(at) * (length(at) - 1) ÷ 2
    pairs = Array{Tuple{Int,Int}, 1}()
    for i in 1:length(at)
        for j in (i+1):length(at)
            push!(pairs, (i, j))
        end
    end
    # @info "num_pairs: $num_pairs, length(pairs): $(length(pairs))"
    energies = Vector{typeof(0.0u"eV")}(undef, length(pairs))
    Threads.@threads for k in eachindex(pairs)
        (i, j) = pairs[k]
        r = pbc_dist(position(at, i), position(at, j), at)
        energies[k] = pair_energy(r, pot)
        # @info "interacting pair: [$(i),$(j)] $(lj_energy(r,lj))"
    end
    return sum(energies)
end

# multi-component interacting energy:
"""
    interacting_energy(at::AbstractSystem, cps::CompositeParameterSets{C,P}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using a multi-component potential (CompositeParameterSets).
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `cps::CompositeParameterSets{C,P}`: The composite potential parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            cps::CompositeParameterSets{C,P}, 
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool}
                            ) where {C,P<:SingleComponentPotential{Pairwise}}
    check_num_components(C, list_num_par, frozen)
    energy = 0.0u"eV"
    components = split_components(at, list_num_par)
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], cps.param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:C
        for j in (i+1):C
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], cps.param_sets[i,j])
            end
        end
    end
    return energy
end


# multi-component interacting energy with an external surface:
"""
    interacting_energy(at::AbstractSystem, cps::CompositeParameterSets{C,P}, list_num_par::Vector{Int}, frozen::Vector{Bool}, surface::AbstractSystem)

Calculate the energy from interactions between free-free and free-frozen particles using a multi-component potential (CompositeParameterSets).
The energy is calculated by summing the pairwise interactions between the free particles and the surface.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `cps::CompositeParameterSets{C,P}`: The composite potential parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.
- `surface::AbstractSystem`: An optional surface system to consider in the energy calculation. See `LJSurfaceWalkers`.  

# Returns
- `energy`: The energy from interactions between particles and the surface.
"""
function interacting_energy(at::AbstractSystem, 
                            cps::CompositeParameterSets{C,P}, 
                            list_num_par::Vector{Int},
                            frozen::Vector{Bool},
                            surface::AbstractSystem
                            ) where {C,P<:SingleComponentPotential{Pairwise}}
    energy = 0.0u"eV"
    components_at = split_components(at, list_num_par)
    components = [components_at..., surface] # combine components from at and surface
    frozen = [frozen..., true] # add frozen state for surface
    list_num_par = [list_num_par..., length(surface)] # combine list_num_par vectors
    num_components = length(components)
    check_num_components(num_components, list_num_par, frozen)
    # @info "num_components: $num_components, list_num_par: $list_num_par, frozen: $frozen"
    # intra-component interactions
    for i in findall(.!frozen) # find non-frozen components
        if length(components[i]) > 1
            energy += intra_component_energy(components[i], cps.param_sets[i,i])
        end
    end
    # inter-component interactions
    for i in 1:num_components
        for j in (i+1):num_components
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], cps.param_sets[i,j])
            end
        end
    end
    return energy
end

# multi-component interacting energy using a single-component potential:
"""
    interacting_energy(at::AbstractSystem, pot::SingleComponentPotential{S}, list_num_par::Vector{Int}, frozen::Vector{Bool})

Calculate the energy from interactions between free-free and free-frozen particles using a single-component potential.
I.e., the components interact with each other using the same parameters.
The energy is calculated by summing the pairwise interactions between the free particles.

# Arguments
- `at::AbstractSystem`: The system for which the energy is calculated.
- `pot::SingleComponentPotential{S}`: The single-component potential.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `frozen::Vector{Bool}`: A vector indicating whether each component is frozen.

# Returns
- `energy`: The energy from interactions between particles.

"""
function interacting_energy(at::AbstractSystem, 
                            pot::SingleComponentPotential{Pairwise},
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
            energy += intra_component_energy(components[i], pot)
        end
    end
    # inter-component interactions
    for i in 1:length(list_num_par)
        for j in (i+1):length(list_num_par)
            if !frozen[i] || !frozen[j] # not both frozen
                # @info "component $i and $j"
                energy += inter_component_energy(components[i], components[j], pot)
            end
        end
    end
    return energy
end
