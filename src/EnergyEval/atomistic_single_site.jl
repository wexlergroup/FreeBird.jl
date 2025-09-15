"""
    single_site_energy(index::Int, at::AbstractSystem, pot::SingleComponentPotential{Pairwise})
    single_site_energy(index::Int, at::AbstractSystem, cps::CompositeParameterSets{C,P}, list_num_par::Vector{Int})
    single_site_energy(index::Int, at::AbstractSystem, cps::CompositeParameterSets{C,P}, list_num_par::Vector{Int}, surface::AbstractSystem)

Calculate the energy of a single site in the system using the Lennard-Jones potential.
The energy is calculated by summing the pairwise interactions between the site and all other sites in the system.

# Arguments
- `index::Int`: The index of the site for which the energy is calculated.
- `at::AbstractSystem`: The system for which the energy is calculated.
- `pot::SingleComponentPotential{Pairwise}`: The Lennard-Jones parameters.
- `cps::CompositeParameterSets{C,P}`: The composite Lennard-Jones parameters.
- `list_num_par::Vector{Int}`: The number of particles in each component.
- `surface::AbstractSystem`: An optional surface system to consider in the energy calculation. See `LJSurfaceWalkers`.

# Returns
- `energy`: The energy of the site.

"""
function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            pot::SingleComponentPotential{Pairwise},
                            list_num_par::Vector{Int}
                            )
    
    all_index = collect(1:length(at))
    popat!(all_index, index)
    energies = Array{typeof(0.0u"eV"), 1}(undef, length(all_index))
    Threads.@threads for i in eachindex(all_index)
        r = pbc_dist(position(at, index), position(at, all_index[i]), at)
        energies[i] = pair_energy(r,pot)
    end
    return sum(energies)
end

function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            cps::CompositeParameterSets{C,P},
                            list_num_par::Vector{Int},
                            ) where {C,P}
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
            energy = pair_energy(r,cps.param_sets[from_comp,from_comp])
            energies[i] = energy
        else
            to_comp = findfirst(x->all_index[i] in x, comp_split)
            energy = pair_energy(r,cps.param_sets[from_comp,to_comp])
            energies[i] = energy
        end
    end
    return sum(energies)
end

function single_site_energy(index::Int,
                            at::AbstractSystem, 
                            cps::CompositeParameterSets{C,P},
                            list_num_par::Vector{Int},
                            surface::AbstractSystem
                            ) where {C,P}
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
            energy = pair_energy(r,cps.param_sets[from_comp,from_comp])
            energies[i] = energy
        else
            to_comp = findfirst(x->all_index[i] in x, comp_split)
            energy = pair_energy(r,cps.param_sets[from_comp,to_comp])
            energies[i] = energy
        end
    end
    internal_e = sum(energies)
    energies_surface = Array{typeof(0.0u"eV"), 1}(undef, length(surface))
    Threads.@threads for i in eachindex(surface.position)
        r = pbc_dist(position(at, index), position(surface, i), at)
        to_comp = C # surface is the last component
        energies_surface[i] = pair_energy(r,cps.param_sets[from_comp,to_comp])
    end
    external_e = sum(energies_surface)
    return internal_e + external_e
end
