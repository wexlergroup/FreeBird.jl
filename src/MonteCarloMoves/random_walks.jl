"""
    single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T

Perform a single atom random walk by updating the position `pos` in each direction by a random amount.
The `step_size` determines the maximum distance the atom can move in any direction.

# Arguments
- `pos::SVector{3,T}`: The current position of the atom as a 3D vector.
- `step_size::Float64`: The maximum distance the atom can move in any direction.

# Returns
- `pos`: The updated position of the atom.

"""
function single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64) where T
    (dx, dy, dz) = (rand(Uniform(-step_size,step_size)) for _ in 1:3) .* unit(T)
    pos = pos .+ (dx, dy, dz)
    return pos 
end

function single_atom_random_walk!(pos::SVector{3,T}, step_size::Float64, dims::Vector{Int}) where T
    ds = [rand(Uniform(-step_size,step_size)) for _ in dims]
    dds = zeros(3)
    for (i, d) in enumerate(dims)
        dds[d] = ds[i]
    end
    ds = dds .* unit(T)
    return pos  .+ ds
end

"""
    single_atom_random_walk!(at::AtomWalker, pot::AbstractPotential, step_size::Float64)
Perform a single atom random walk on the `AtomWalker` object `at` using the potential `pot` and a specified `step_size`.

# Arguments
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `pot::AbstractPotential`: The potential energy function for the system.
- `step_size::Float64`: The maximum distance an atom can move in any direction.

# Returns
- `at::AtomWalker`: The updated walker after the random walk.
"""
function single_atom_random_walk!(at::AtomWalker{C}, pot::AbstractPotential, step_size::Float64) where C
    config = at.configuration
    # select a random free atom to move
    free_index = free_par_index(at)
    i_at = rand(free_index)
    # calculate the energy before the move
    prewalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
    # prewalk_potential = interacting_energy(config, lj, at.list_num_par, at.frozen)
    # get the current position of the atom
    pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
    # perform the random walk
    pos = single_atom_random_walk!(pos, step_size)
    # wrap the atom around the periodic boundary
    pos = periodic_boundary_wrap!(pos, config)
    # update the position of the atom
    config.position[i_at] = pos
    # calculate the energy after the move
    postwalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
    # postwalk_potential = interacting_energy(config, lj, at.list_num_par, at.frozen)
    # calculate the energy difference
    e_diff = postwalk_energy - prewalk_energy
    # e_pot_diff = postwalk_potential - prewalk_potential

    return at, e_diff
end

"""
    MC_random_walk!(n_steps::Int, at::AtomWalker, pot::AbstractPotential, step_size::Float64, emax::typeof(0.0u"eV"))

Perform a Monte Carlo random walk on the atomic/molecular system.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `pot::AbstractPotential`: The potential energy function for the system.
- `step_size::Float64`: The maximum distance an atom can move in any direction.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `at::AtomWalker`: The updated walker.

"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::AbstractPotential, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV")
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        # prewalk_energy = interacting_energy(config, pot, at.list_num_par, at.frozen)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        postwalk_energy = interacting_energy(config, pot, at.list_num_par, at.frozen)

        if postwalk_energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = postwalk_energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end

function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::SingleComponentPotential{ManyBody}, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV")
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        # prewalk_energy = interacting_energy(config, pot, at.list_num_par, at.frozen)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        postwalk_energy = interacting_energy(config, pot)

        if postwalk_energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = postwalk_energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end

"""
    MC_random_walk!(n_steps::Int, at::AtomWalker, pot::LennardJonesParameterSets, step_size::Float64, emax::typeof(0.0u"eV"))

Perform a Monte Carlo random walk on the atomic/molecular system. Specialized for Lennard-Jones potentials.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `pot::LennardJonesParameterSets`: The potential energy function for the system.
- `step_size::Float64`: The maximum distance an atom can move in any direction.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `at::AtomWalker`: The updated walker.

"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::LennardJonesParameterSets, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV")
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        prewalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        postwalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
        e_diff = postwalk_energy - prewalk_energy
        # energy = interacting_energy(config, lj, at.list_num_par, at.frozen) + at.energy_frozen_part
        energy = at.energy + e_diff
        if energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end

"""    MC_random_walk!(n_steps::Int, at::AtomWalker{C}, pot::AbstractPotential, step_size::Float64, emax::typeof(0.0u"eV"), surface::AtomWalker{CS})
Perform a Monte Carlo random walk on the atomic/molecular system with an external surface.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `pot::AbstractPotential`: The potential energy function for the system.
- `step_size::Float64`: The maximum distance an atom can move in any direction.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.
- `surface::AtomWalker{CS}`: The surface walker object to consider in the energy calculation. Typically frozen. 

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
"""
function MC_random_walk!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::AbstractPotential, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV"),
                    surface::AtomWalker{CS},
                    ) where {C, CS}
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        prewalk_energy = single_site_energy(i_at, config, pot, at.list_num_par, surface.configuration)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        postwalk_energy = single_site_energy(i_at, config, pot, at.list_num_par, surface.configuration)
        e_diff = postwalk_energy - prewalk_energy
        # energy = interacting_energy(config, lj, at.list_num_par, at.frozen) + at.energy_frozen_part
        energy = at.energy + e_diff
        if energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end

"""
    MC_random_walk_2D!(n_steps::Int, at::AtomWalker, pot::AbstractPotential, step_size::Float64, emax::typeof(0.0u"eV"); dims::Vector{Int}=[1,2])

Perform a Monte Carlo random walk on the atomic/molecular system in 2D.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `pot::AbstractPotential`: The potential energy function for the system.
- `step_size::Float64`: The maximum distance an atom can move in any direction.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.
- `dims::Vector{Int}=[1,2]`: The dimensions in which the random walk is performed.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `at::AtomWalker`: The updated walker.

"""
function MC_random_walk_2D!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::AbstractPotential, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV");
                    dims::Vector{Int}=[1,2]
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size, dims)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        energy = interacting_energy(config, pot, at.list_num_par, at.frozen) + at.energy_frozen_part
        if energy >= emax
            # reject the move, revert to original position
            config.position[i_at] = orig_pos
        else
            at.energy = energy
            # accept the move
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end

"""
    MC_random_walk!(n_steps::Int, lattice::LatticeWalker, h::ClassicalHamiltonian, emax::Float64; energy_perturb::Float64=0.0)

Perform a Monte Carlo random walk on the lattice system.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `lattice::LatticeWalker`: The walker to perform the random walk on.
- `h::ClassicalHamiltonian`: The lattice gas Hamiltonian.
- `emax::Float64`: The maximum energy allowed for accepting a move.
- `energy_perturb::Float64=0.0`: The energy perturbation used to make degenerate configurations distinguishable.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `lattice::LatticeWalker`: The updated walker.

"""
function MC_random_walk!(n_steps::Int,
                         lattice::LatticeWalker{C},
                         h::ClassicalHamiltonian,
                         emax::Float64;
                         energy_perturb::Float64=0.0,
                         ) where C

    n_accept = 0
    accept_this_walker = false
    emax = emax * unit(lattice.energy)

    for i_mc_step in 1:n_steps
        current_lattice = lattice.configuration

        proposed_lattice = deepcopy(current_lattice)

        lattice_random_walk!(proposed_lattice)
        
        perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
        proposed_energy = interacting_energy(proposed_lattice, h) + perturbation_energy

        @debug "proposed_energy = $proposed_energy, perturbed_energy = $(perturbation_energy), emax = $(emax)), accept = $(proposed_energy < emax)"
        if proposed_energy >= emax
            continue
        else
            lattice.configuration = proposed_lattice
            lattice.energy = proposed_energy
            n_accept += 1
            accept_this_walker = true
        end
    end
    return accept_this_walker, n_accept/n_steps, lattice
end

"""
    MC_new_sample!(lattice::LatticeWalker, h::ClassicalHamiltonian, emax::Float64; energy_perturb::Float64=0.0)

Generate a new sample for the lattice system.

# Arguments
- `lattice::LatticeWalker`: The walker to generate a new sample for.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `emax::Float64`: The maximum energy allowed for accepting a move.
- `energy_perturb::Float64=0.0`: The energy perturbation used to make degenerate configurations distinguishable.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `lattice::LatticeWalker`: The updated walker.

"""
function MC_new_sample!(lattice::LatticeWalker{C},
                        h::ClassicalHamiltonian,
                        emax::Float64;
                        energy_perturb::Float64=0.0,
                        ) where C

    accept_this_walker = false
    emax = emax * unit(lattice.energy)

    current_lattice = lattice.configuration
    proposed_lattice = deepcopy(current_lattice)
    generate_random_new_lattice_sample!(proposed_lattice)

    perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
    proposed_energy = interacting_energy(proposed_lattice, h) + perturbation_energy

    @debug "proposed_energy = $proposed_energy, perturbed_energy = $(perturbation_energy), emax = $(emax)), accept = $(proposed_energy < emax)"

    if proposed_energy < emax
        lattice.configuration = proposed_lattice
        lattice.energy = proposed_energy
        accept_this_walker = true
    end

    return accept_this_walker, lattice
end

"""
    MC_rejection_sampling!(lattice::LatticeWalker, h::ClassicalHamiltonian, emax::Float64; energy_perturb::Float64=0.0, max_iter=10_000)

Perform a Monte Carlo rejection sampling on the lattice system.

# Arguments
- `lattice::LatticeWalker`: The walker to perform the rejection sampling on.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `emax::Float64`: The maximum energy allowed for accepting a move.
- `energy_perturb::Float64=0.0`: The energy perturbation used to make degenerate configurations distinguishable.
- `max_iter::Int=10_000`: The maximum number of iterations to perform.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `lattice::LatticeWalker`: The updated walker.

"""
function MC_rejection_sampling!(lattice::LatticeWalker{C},
                        h::ClassicalHamiltonian,
                        emax::Float64;
                        energy_perturb::Float64=0.0,
                        max_iter::Int = 10_000,
                        ) where C

    accept_this_walker = false
    emax = emax * unit(lattice.energy)

    current_energy = lattice.energy # initialize to a value greater than emax
    current_lattice = lattice.configuration

    counter = 0

    while current_energy >= emax
        
        proposed_lattice = deepcopy(current_lattice)
        generate_random_new_lattice_sample!(proposed_lattice)

        perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
        raw_energy = interacting_energy(proposed_lattice, h)
        proposed_energy = raw_energy + perturbation_energy

        counter += 1
        if counter > max_iter
            accept_this_walker = false
            break
        end

        @debug "proposed_energy = $proposed_energy, perturbed_energy = $(perturbation_energy), emax = $(emax)), accept = $(proposed_energy < emax)"

        if proposed_energy <= emax
            lattice.configuration = proposed_lattice
            lattice.energy = proposed_energy
            accept_this_walker = true
            break
        end
    end

    return accept_this_walker, lattice
end


"""
    generate_random_new_lattice_sample!(lattice::MLattice{C}) where C

Generate a new random sample for the multi-component lattice system.

# Arguments
- `lattice::MLattice{C}`: The lattice system to generate a new sample for.

# Returns
- `lattice::MLattice{C}`: The updated lattice system.

"""
function generate_random_new_lattice_sample!(lattice::MLattice{C}) where C
    occupancy = occupied_site_count(lattice)
    # flush occupancy
    for i in eachindex(lattice.components)
        for j in eachindex(lattice.components[i])
            lattice.components[i][j] = false
        end
    end
    # generate a "look-up" table of unoccupied sites with elements get deleted as sites are occupied
    # faster than shuffling a list of components
    unoccupied = collect(1:num_sites(lattice))
    for i in eachindex(lattice.components)
        samples = sample(unoccupied, occupancy[i], replace=false, ordered=true)
        for j in samples
            lattice.components[i][j] = true
            deleteat!(unoccupied, findfirst(isequal(j), unoccupied))
        end
    end
    return lattice
end

"""
    generate_random_new_lattice_sample!(lattice::SLattice)

Generate a new random sample for the single-component lattice system.
"""
function generate_random_new_lattice_sample!(lattice::SLattice)
    number_occupied_sites = sum(lattice.components[1])
    # flush occupancy
    lattice.components[1] .= false
    for i in sample(eachindex(lattice.components[1]), number_occupied_sites, replace=false)
        lattice.components[1][i] = true
    end
    return lattice
end

"""
    lattice_random_walk!(lattice::SLattice)  

Perform a Monte Carlo random walk on the single-component lattice system.

# Arguments
- `lattice::SLattice`: The single-component lattice system to perform the random walk on.
# Returns
- `lattice::SLattice`: The proposed lattice after the random walk.
"""
function lattice_random_walk!(lattice::SLattice)
    # pick a random site to hop from
    hop_from = rand(eachindex(lattice.components[1]))
    # pick a random site to hop to (can be the same as hop_from)
    hop_to = rand(eachindex(lattice.components[1]))
    # propose a swap in occupation state (only if it maintains constant N)
    # proposed_lattice = deepcopy(lattice)
    if lattice.components[1][hop_from] != lattice.components[1][hop_to]
        lattice.components[1][hop_from], lattice.components[1][hop_to] = 
        lattice.components[1][hop_to], lattice.components[1][hop_from]
    end
    return lattice
end

"""
    lattice_random_walk!(lattice::MLattice{C,G}) where {C,G}

Perform a Monte Carlo random walk on the multi-component lattice system.

# Arguments
- `lattice::MLattice{C,G}`: The multi-component lattice system to perform the random walk on.
# Returns
- `lattice::MLattice{C,G}`: The proposed lattice after the random walk.
"""
function lattice_random_walk!(lattice::MLattice{C,G}) where {C,G}
    # pick a random component to hop in
    picked_comp::Int = rand(1:C)
    # pick a random site to hop from
    hop_from::Int = rand(eachindex(lattice.components[picked_comp]))
    # pick a random site to hop to (can be the same as hop_from)
    hop_to::Int = rand(eachindex(lattice.components[picked_comp]))

    # println("hop from: $hop_from, hop to: $hop_to, comp: $comp") # debug
    
    # case 1: occupied/unoccupied swap
    is_occupied_from::Bool = any([lattice.components[comp][hop_from] for comp in 1:C])
    is_occupied_to::Bool = any([lattice.components[comp][hop_to] for comp in 1:C])
    # println("is_occupied_from: $is_occupied_from, is_occupied_to: $is_occupied_to") # debug
    if is_occupied_from != is_occupied_to # only swap if the occupation state changes
        # swap the occupation state of the sites
        return swap_empty_occupied_sites!(lattice, hop_from, hop_to)
    # case 2: both sites occupied, swap components
    elseif is_occupied_from && is_occupied_to
        return swap_occupied_sites_across_components!(lattice, hop_from, hop_to)
    end
    # case 3: both sites unoccupied, do nothing
    return lattice
end

"""
    swap_empty_occupied_sites!(lattice::MLattice{C,G}, hop_from::Int, hop_to::Int) where {C,G}

Swap the occupation state of two sites in the lattice of any component.

# Arguments
- `lattice::MLattice{C,G}`: The lattice to perform the swap on.
- `hop_from::Int`: The index of the site to hop from.
- `hop_to::Int`: The index of the site to hop to.
# Returns
- `lattice::MLattice{C,G}`: The updated lattice after the swap.
"""
function swap_empty_occupied_sites!(lattice::MLattice{C,G}, 
                                     hop_from::Int, 
                                     hop_to::Int) where {C,G}
    for comp::Int in 1:C
        lattice.components[comp][hop_from], lattice.components[comp][hop_to] = 
        lattice.components[comp][hop_to], lattice.components[comp][hop_from]
    end
    return lattice
end

"""    swap_occupied_sites_across_components!(lattice::MLattice{C,G}, hop_from::Int, hop_to::Int) where {C,G}
Swap the occupation state of two sites across different components in the lattice.
# Arguments
- `lattice::MLattice{C,G}`: The lattice to perform the swap on.
- `hop_from::Int`: The index of the site to hop from.
- `hop_to::Int`: The index of the site to hop to.
# Returns
- `lattice::MLattice{C,G}`: The updated lattice after the swap.
"""
function swap_occupied_sites_across_components!(lattice::MLattice{C,G}, 
                                     hop_from::Int, 
                                     hop_to::Int) where {C,G}
    # find out which components are occupied at the sites
    comp_from::Int = findfirst(isequal(true), [lattice.components[comp][hop_from] for comp in 1:C])
    comp_to::Int = findfirst(isequal(true), [lattice.components[comp][hop_to] for comp in 1:C])
    if comp_from != comp_to # only swap if the components are different
        # first, unoccupy the site at comp_from
        lattice.components[comp_from][hop_from] = !lattice.components[comp_from][hop_from]
        # second, unoccupy the site at comp_to
        lattice.components[comp_to][hop_from] = !lattice.components[comp_to][hop_from]
        # third, occupy the site at comp_to
        lattice.components[comp_to][hop_to] = !lattice.components[comp_to][hop_to]
        # finally, occupy the site at comp_from
        lattice.components[comp_from][hop_to] = !lattice.components[comp_from][hop_to]
    end
    return lattice
end