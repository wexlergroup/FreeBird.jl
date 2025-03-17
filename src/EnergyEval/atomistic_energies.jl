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
    pbc = periodicity(at)
    box = cell_vectors(at)
    distsq = 0.0u"Å"^2
    for i in eachindex(pos1)
        if pbc[i] == true
            distsq += min(abs(pos1[i] - pos2[i]), box[i][i] - abs(pos1[i] - pos2[i]))^2
        elseif pbc[i] == false
            distsq += (pos1[i] - pos2[i])^2
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
        energies[k] = lj_energy(r,lj)
        # @info "interacting pair: [$(i),$(j)] $(lj_energy(r,lj))"
    end
    return sum(energies)
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
    single_site_energy(index::Int, at::AbstractSystem, lj::LJParameters)
    single_site_energy(index::Int, at::AbstractSystem, ljs::CompositeLJParameters{C}, list_num_par::Vector{Int})

Calculate the energy of a single site in the system using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the site and all other sites in the system.

# Arguments
- `index::Int`: The index of the site for which the energy is calculated.
- `at::AbstractSystem`: The system for which the energy is calculated.
- `lj::LJParameters`: The Lennard-Jones parameters.
- `ljs::CompositeLJParameters{C}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.

# Returns
- `energy`: The energy of the site.

"""
function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            lj::LJParameters,
                            list_num_par::Vector{Int}
                            )
    
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        energies[i] = lj_energy(r,lj)
    end
    return sum(energies)
end

function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            ljs::CompositeLJParameters{C},
                            list_num_par::Vector{Int},
                            ) where {C}
    comp_cut = vcat([0],cumsum(list_num_par))
    # @info "comp_cut: $comp_cut"
    comp_split = [comp_cut[i]+1:comp_cut[i+1] for i in 1:length(list_num_par)]
    # @info "comp_split: $comp_split"
    from_comp = findfirst(x->index in x, comp_split)
    # @info "from_comp: $from_comp"
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        if all_index[i] in comp_split[from_comp]
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,from_comp])
            energies[i] = energy
        else
            to_comp = findfirst(x->all_index[i] in x, comp_split)
            energy = lj_energy(r,ljs.lj_param_sets[from_comp,to_comp])
            energies[i] = energy
        end
    end
    return sum(energies)
end
