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
        config = lattice.configuration

        # In-place proposal: the hop-pair walk mutates the live configuration
        # and returns the drawn pair; a rejected proposal is reverted by
        # replaying the pair (involution). No random draw changes.
        hop_from, hop_to = _lattice_walk_draw!(config)

        perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
        proposed_energy = interacting_energy(config, h) + perturbation_energy

        @debug "proposed_energy = $proposed_energy, perturbed_energy = $(perturbation_energy), emax = $(emax)), accept = $(proposed_energy < emax)"
        if proposed_energy >= emax
            _lattice_walk_apply!(config, hop_from, hop_to)
            continue
        else
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
    _lattice_walk_draw!(lattice)
    return lattice
end

"""
    _lattice_walk_draw!(lattice) -> (hop_from, hop_to)

Internal: perform one hop-pair walk step with exactly the draws of
`lattice_random_walk!`, in the same order, returning the drawn pair so the
copy-free kernels can revert a rejected proposal. The conditional exchange
applied for a drawn pair is an involution, so replaying the pair through
`_lattice_walk_apply!` reverts the step exactly.
"""
function _lattice_walk_draw!(lattice::SLattice)
    # pick a random site to hop from
    hop_from = rand(eachindex(lattice.components[1]))
    # pick a random site to hop to (can be the same as hop_from)
    hop_to = rand(eachindex(lattice.components[1]))
    _lattice_walk_apply!(lattice, hop_from, hop_to)
    return hop_from, hop_to
end

function _lattice_walk_draw!(lattice::MLattice{C,G}) where {C,G}
    # pick a random component to hop in (draw kept for stream identity;
    # all component vectors share the same site index range)
    picked_comp::Int = rand(1:C)
    # pick a random site to hop from
    hop_from::Int = rand(eachindex(lattice.components[picked_comp]))
    # pick a random site to hop to (can be the same as hop_from)
    hop_to::Int = rand(eachindex(lattice.components[picked_comp]))
    _lattice_walk_apply!(lattice, hop_from, hop_to)
    return hop_from, hop_to
end

"""
    _lattice_walk_apply!(lattice, hop_from, hop_to)

Internal: apply the conditional occupancy exchange of the hop-pair walk for an
already-drawn pair, with no random draws. Self-inverse for every case (the
empty/occupied exchange, the cross-component exchange, and the no-op), so a
second call with the same pair reverts the first.
"""
function _lattice_walk_apply!(lattice::SLattice, hop_from::Int, hop_to::Int)
    # propose a swap in occupation state (only if it maintains constant N)
    if lattice.components[1][hop_from] != lattice.components[1][hop_to]
        lattice.components[1][hop_from], lattice.components[1][hop_to] =
        lattice.components[1][hop_to], lattice.components[1][hop_from]
    end
    return lattice
end

function _lattice_walk_apply!(lattice::MLattice{C,G}, hop_from::Int, hop_to::Int) where {C,G}
    # case 1: occupied/unoccupied swap
    is_occupied_from::Bool = any([lattice.components[comp][hop_from] for comp in 1:C])
    is_occupied_to::Bool = any([lattice.components[comp][hop_to] for comp in 1:C])
    if is_occupied_from != is_occupied_to # only swap if the occupation state changes
        return swap_empty_occupied_sites!(lattice, hop_from, hop_to)
    # case 2: both sites occupied, swap components
    elseif is_occupied_from && is_occupied_to
        return swap_occupied_sites_across_components!(lattice, hop_from, hop_to)
    end
    # case 3: both sites unoccupied, do nothing
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
    _lattice_walk_draw!(lattice)
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


"""
    MC_cluster_walk!(n_steps::Int, lattice::LatticeWalker{C}, h::ClassicalHamiltonian, emax::Float64, cluster_p::Float64; energy_perturb::Float64=0.0)

Perform a sequence of geometric cluster moves on the lattice system, accepting each if `E < emax`.

# Arguments
- `n_steps::Int`: The number of cluster move attempts to perform.
- `lattice::LatticeWalker{C}`: The walker to perform cluster moves on.
- `h::ClassicalHamiltonian`: The lattice Hamiltonian.
- `emax::Float64`: The maximum energy allowed for accepting a move (dimensionless).
- `cluster_p::Float64`: The growth probability for BFS cluster construction.
- `energy_perturb::Float64=0.0`: Energy perturbation to break degeneracies.

# Returns
- `accept_this_walker::Bool`: Whether at least one move was accepted.
- `accept_rate::Float64`: The fraction of accepted moves.
- `lattice::LatticeWalker`: The updated walker.
"""
function MC_cluster_walk!(n_steps::Int,
                          lattice::LatticeWalker{C},
                          h::ClassicalHamiltonian,
                          emax::Float64,
                          cluster_p::Float64;
                          energy_perturb::Float64=0.0) where C
    n_accept = 0
    accept_this_walker = false
    emax_u = emax * unit(lattice.energy)

    cluster_pairs = Tuple{Int,Int}[]
    for _ in 1:n_steps
        config = lattice.configuration
        # In-place proposal: record the applied pairs so a rejected cluster
        # move reverts by re-applying them (involution). No random draw changes.
        empty!(cluster_pairs)
        geometric_cluster_swap!(config, cluster_p; record=cluster_pairs)

        perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
        proposed_energy = interacting_energy(config, h) + perturbation_energy

        if proposed_energy < emax_u
            lattice.energy = proposed_energy
            n_accept += 1
            accept_this_walker = true
        else
            _apply_cluster_pairs!(config, cluster_pairs)
        end
    end
    return accept_this_walker, n_accept / max(n_steps, 1), lattice
end


# ======================================================================
# Grand-canonical move primitives
# ======================================================================

"""
    random_microstate!(lattice::SLattice; p::Float64=0.5)

Set each site occupied independently with probability `p`, producing a
variable-N configuration suitable for grand-canonical sampling.

# Arguments
- `lattice::SLattice`: The single-component lattice to randomize.
- `p::Float64=0.5`: Per-site occupation probability.

# Returns
- `lattice::SLattice`: The mutated lattice with a random microstate.
"""
function random_microstate!(lattice::SLattice; p::Float64=0.5)
    for i in eachindex(lattice.components[1])
        lattice.components[1][i] = rand() < p
    end
    return lattice
end

"""
    lattice_insert_particle!(lattice::SLattice)

Insert a particle at a random empty site. Returns `true` if successful,
`false` if the lattice is full.

# Arguments
- `lattice::SLattice`: The single-component lattice.

# Returns
- `success::Bool`: Whether a particle was inserted.
- `lattice::SLattice`: The mutated lattice.
- `site::Int`: The inserted site index (0 when unsuccessful). A trailing
  addition: two-name destructures of the previous return keep working.
"""
function lattice_insert_particle!(lattice::SLattice)
    n_sites = num_sites(lattice)
    n_occ = sum(lattice.components[1])
    if n_occ >= n_sites
        return false, lattice, 0
    end
    # Collect empty site indices
    empty_sites = findall(.!lattice.components[1])
    site = rand(empty_sites)
    lattice.components[1][site] = true
    return true, lattice, site
end

"""
    lattice_delete_particle!(lattice::SLattice)

Delete a particle from a random occupied site. Returns `true` if successful,
`false` if the lattice is empty.

# Arguments
- `lattice::SLattice`: The single-component lattice.

# Returns
- `success::Bool`: Whether a particle was deleted.
- `lattice::SLattice`: The mutated lattice.
- `site::Int`: The vacated site index (0 when unsuccessful). A trailing
  addition: two-name destructures of the previous return keep working.
"""
function lattice_delete_particle!(lattice::SLattice)
    n_occ = sum(lattice.components[1])
    if n_occ == 0
        return false, lattice, 0
    end
    occupied_sites = findall(lattice.components[1])
    site = rand(occupied_sites)
    lattice.components[1][site] = false
    return true, lattice, site
end

"""
    lattice_biased_sites(lattice::SLattice; predicate::Symbol=:contact, shells::Int=1)

Return the indices of empty sites selected by an occupancy predicate over
neighbor shells `1:shells`:

- `:contact`: empty sites with at least one occupied neighbor.
- `:cavity`: empty sites with no occupied neighbor.

For fixed `shells` the two predicates partition the empty sites:
`length(contact set) + length(cavity set) == M - N`. Neighbor lists carrying
periodic-image multiplicity (the same neighbor listed more than once) work
as-is: any occupied appearance marks contact.

Useful as a nested-sampling observable, e.g.
`:n_cavity => cfg -> length(lattice_biased_sites(cfg; predicate=:cavity))`.

Throws `ArgumentError` for an unknown predicate, `shells < 1`, or `shells`
exceeding the lattice's neighbor-shell count.
"""
function lattice_biased_sites(lattice::SLattice; predicate::Symbol=:contact, shells::Int=1)
    if predicate !== :contact && predicate !== :cavity
        throw(ArgumentError("unknown predicate :$predicate; expected :contact or :cavity"))
    end
    if shells < 1
        throw(ArgumentError("shells must be >= 1, got $shells"))
    end
    neighbors = lattice.neighbors
    n_shells = isempty(neighbors) ? 0 : length(neighbors[1])
    if shells > n_shells
        throw(ArgumentError(
            "the biased-site predicate spans $shells neighbor shells but the " *
            "lattice provides only $n_shells (= length(cutoff_radii)); " *
            "extend cutoff_radii so every counted shell exists"))
    end
    occ = lattice.components[1]
    sites = Int[]
    for site in eachindex(occ)
        occ[site] && continue
        has_occupied_neighbor = false
        for shell in 1:shells
            for nb in neighbors[site][shell]
                if occ[nb]
                    has_occupied_neighbor = true
                    break
                end
            end
            has_occupied_neighbor && break
        end
        if predicate === :contact ? has_occupied_neighbor : !has_occupied_neighbor
            push!(sites, site)
        end
    end
    return sites
end

"""
    MC_grand_canonical_walk!(n_steps::Int, lattice::LatticeWalker{1},
                             h::ClassicalHamiltonian, omega_max::Float64,
                             mu::Float64;
                             p_move::Float64=0.5, p_insert::Float64=0.25,
                             energy_perturb::Float64=0.0, n_max::Int=typemax(Int),
                             clusters_freq::Int=0, swaps_freq::Int=1,
                             cluster_p::Float64=0.3, z0::Float64=1.0,
                             p_bias::Float64=0.0, bias_predicate::Symbol=:contact,
                             bias_shells::Int=1)

Perform grand-canonical MCMC on a single-component lattice, mixing fixed-N
moves (local swaps and/or geometric cluster moves) with single-site insertion
and deletion.

Each step:
- With probability `p_move`: propose a fixed-N move (local swap or cluster move,
  selected by `swaps_freq:clusters_freq` ratio).
- With probability `p_insert`: propose inserting one particle.
- With probability `1 - p_move - p_insert`: propose deleting one particle.

Insert/delete proposals use a Metropolis correction to preserve the ideal-lattice-gas
prior at reference fugacity `z0` — the Bernoulli product measure giving each
configuration with N particles weight `z0^N` — over all microstates with Ω < Ω_max:
- Insert ratio: `z0 * (p_delete / p_insert) * (M - N) / (N + 1)`
- Delete ratio: `(1 / z0) * (p_insert / p_delete) * N / (M - N + 1)`

The default `z0 = 1.0` reduces to the uniform prior over all 2^M microstates,
matching the Ω-sorted grand-canonical nested sampling construction. `z0 ≠ 1`
is used by the ideal-gas-referenced (E-sorted) construction, where the walk
runs with `mu = 0` so the Ω ceiling reduces to an energy ceiling.

Cluster moves are symmetric (no Metropolis correction), accepted if Ω < Ω_max.

# Arguments
- `n_steps::Int`: Number of MCMC steps.
- `lattice::LatticeWalker{1}`: The walker (single component).
- `h::ClassicalHamiltonian`: The lattice Hamiltonian.
- `omega_max::Float64`: Upper bound on grand potential Ω = E − μN (unitless).
- `mu::Float64`: Chemical potential (unitless, in same energy units as Hamiltonian).
- `p_move::Float64=0.5`: Probability of a fixed-N move.
- `p_insert::Float64=0.25`: Probability of an insertion move.
- `energy_perturb::Float64=0.0`: Energy perturbation for degeneracy breaking.
- `n_max::Int=typemax(Int)`: Upper bound on particle count.
- `clusters_freq::Int=0`: Relative weight of cluster moves within fixed-N branch (0 = disabled).
- `swaps_freq::Int=1`: Relative weight of local swaps within fixed-N branch.
- `cluster_p::Float64=0.3`: Current cluster growth probability.
- `z0::Float64=1.0`: Reference fugacity of the prior preserved by insert/delete
  (1.0 = uniform prior over microstates).
- `p_bias::Float64=0.0`: Probability that an insertion draws its site uniformly
  from `lattice_biased_sites(x; predicate=bias_predicate, shells=bias_shells)`
  instead of from all empty sites. The acceptance then uses the composite
  proposal density evaluated for the actually chosen site, and the deletion
  acceptance uses the reverse composite density evaluated on the post-deletion
  configuration (every set quantity on the lower-N member of the pair; the
  delete ratio's particle count `n` is the pre-delete count), so the
  `z0`-weighted prior stays invariant. An empty biased set makes the biased
  sub-channel a null proposal (counted as attempted, nothing changes); a zero
  reverse density (`p_bias = 1` with the vacated site outside the biased set)
  is an immediate reject. `0.0` reproduces the legacy sampler bit-for-bit.
- `bias_predicate::Symbol=:contact`: Biased-set predicate, `:contact` or `:cavity`.
- `bias_shells::Int=1`: Neighbor shells scanned by the predicate.

# Returns
- `accept_this_walker::Bool`: Whether at least one move was accepted.
- `accept_rate::Float64`: Fraction of accepted moves.
- `lattice::LatticeWalker{1}`: The updated walker.
- `cluster_accepted_count::Int`: Number of accepted cluster moves (for adaptive tuning).
- `cluster_total_count::Int`: Number of attempted cluster moves (for adaptive tuning).
- `move_stats::NamedTuple`: Per-move-type attempt/accept counters for the walk:
  `swap_*`, `cluster_*` (duplicating the two preceding elements),
  `insert_uniform_*`, `insert_biased_*`, `delete_*`. Attempts are counted at
  proposal: ceiling rejections included, guard skips excluded, and an
  empty-biased-set null proposal counts as an attempted biased insert.
"""
function MC_grand_canonical_walk!(n_steps::Int,
                                  lattice::LatticeWalker{1},
                                  h::ClassicalHamiltonian,
                                  omega_max::Float64,
                                  mu::Float64;
                                  p_move::Float64=0.5,
                                  p_insert::Float64=0.25,
                                  energy_perturb::Float64=0.0,
                                  n_max::Int=typemax(Int),
                                  clusters_freq::Int=0,
                                  swaps_freq::Int=1,
                                  cluster_p::Float64=0.3,
                                  z0::Float64=1.0,
                                  p_bias::Float64=0.0,
                                  bias_predicate::Symbol=:contact,
                                  bias_shells::Int=1)
    if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
        throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
    end
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive"))
    end
    if !(0.0 <= p_bias <= 1.0)
        throw(ArgumentError("p_bias must satisfy 0 <= p_bias <= 1"))
    end
    if bias_predicate !== :contact && bias_predicate !== :cavity
        throw(ArgumentError("unknown bias_predicate :$bias_predicate; expected :contact or :cavity"))
    end
    if bias_shells < 1
        throw(ArgumentError("bias_shells must be >= 1, got $bias_shells"))
    end
    if p_bias > 0.0
        nbrs = lattice.configuration.neighbors
        n_lattice_shells = isempty(nbrs) ? 0 : length(nbrs[1])
        if bias_shells > n_lattice_shells
            throw(ArgumentError(
                "the biased insertion channel spans $bias_shells neighbor shells but the " *
                "lattice provides only $n_lattice_shells (= length(cutoff_radii)); " *
                "extend cutoff_radii so every counted shell exists"))
        end
    end

    n_accept = 0
    accept_this_walker = false
    p_delete = 1.0 - p_move - p_insert
    n_sites = num_sites(lattice.configuration)
    n_cap = min(n_sites, n_max)
    omega_max_u = omega_max * unit(lattice.energy)

    # Cluster move mixing: compute probability of cluster vs local swap within fixed-N branch
    total_fixed_n_freq = clusters_freq + swaps_freq
    p_cluster_in_move = total_fixed_n_freq > 0 ? clusters_freq / total_fixed_n_freq : 0.0

    cluster_accepted_count = 0
    cluster_total_count = 0
    swap_attempted = 0
    swap_accepted = 0
    insert_uniform_attempted = 0
    insert_uniform_accepted = 0
    insert_biased_attempted = 0
    insert_biased_accepted = 0
    delete_attempted = 0
    delete_accepted = 0

    # Reused inside the loop only when p_bias > 0 (the default path never
    # computes biased sets)
    biased_set = Int[]
    insert_site = 0
    insert_from_biased = false
    deleted_site = 0
    # In-place proposal state: proposals mutate the live configuration and a
    # rejected proposal is reverted from this per-move record (hop pair,
    # cluster pair list, or flipped site). Every revert is an involution
    # replay or a single flip, so the reverted state is bit-identical, and
    # no random draw is added or removed anywhere on this path.
    hop_from = 0
    hop_to = 0
    cluster_pairs = Tuple{Int,Int}[]

    for _ in 1:n_steps
        r = rand()
        config = lattice.configuration
        n = sum(config.components[1])

        if r < p_move
            # Fixed-N branch: choose cluster or local swap
            if p_cluster_in_move > 0.0 && rand() < p_cluster_in_move
                # Geometric cluster move
                empty!(cluster_pairs)
                geometric_cluster_swap!(config, cluster_p; record=cluster_pairs)
                move_type = :cluster
                cluster_total_count += 1
            else
                # Local swap
                hop_from, hop_to = _lattice_walk_draw!(config)
                move_type = :move
                swap_attempted += 1
            end
        elseif r < p_move + p_insert
            # Insertion
            if n >= n_cap || p_insert <= 0.0
                continue                # guard skip: not counted as an attempt
            end
            if p_bias > 0.0
                # Composite channel. S(x) is computed ONCE per attempt and
                # reused for both the biased draw and the composite proposal
                # density in the acceptance below.
                biased_set = lattice_biased_sites(config;
                                                  predicate=bias_predicate,
                                                  shells=bias_shells)
                if rand() < p_bias
                    insert_biased_attempted += 1
                    if isempty(biased_set)
                        continue        # null proposal: counted, no further draws
                    end
                    insert_site = rand(biased_set)
                    insert_from_biased = true
                else
                    # Uniform sub-channel: the same single rand(::Vector) draw
                    # as lattice_insert_particle!, inlined to capture the site
                    # for the composite density (keep in lockstep with it)
                    insert_uniform_attempted += 1
                    empty_sites = findall(.!config.components[1])
                    insert_site = rand(empty_sites)
                    insert_from_biased = false
                end
                config.components[1][insert_site] = true
            else
                # p_bias == 0: legacy path, bit-identical RNG stream (no
                # channel draw; lattice_insert_particle! draws exactly once)
                success, _, insert_site = lattice_insert_particle!(config)
                if !success
                    continue
                end
                insert_uniform_attempted += 1
                insert_from_biased = false
            end
            move_type = :insert
        else
            # Deletion: site selection is uniform over occupied sites in both
            # paths; only the acceptance differs
            if n == 0 || p_delete <= 0.0
                continue                # guard skip: not counted as an attempt
            end
            if p_bias > 0.0
                # Inlined uniform deletion (the same single rand(::Vector)
                # draw as lattice_delete_particle!) to capture the vacated
                # site for the reverse composite density (keep in lockstep)
                occupied_sites = findall(config.components[1])
                deleted_site = rand(occupied_sites)
                config.components[1][deleted_site] = false
            else
                success, _, deleted_site = lattice_delete_particle!(config)
                if !success
                    continue
                end
            end
            delete_attempted += 1
            move_type = :delete
        end

        perturbation_energy = energy_perturb * (rand() - 0.5) * unit(lattice.energy)
        proposed_energy = interacting_energy(config, h) + perturbation_energy
        n_new = sum(config.components[1])
        proposed_omega = proposed_energy - mu * n_new * unit(lattice.energy)

        if proposed_omega >= omega_max_u
            _gc_revert_move!(config, move_type, hop_from, hop_to,
                             cluster_pairs, insert_site, deleted_site)
            continue
        end

        # Metropolis correction for insert/delete detailed balance under the
        # z0^N-weighted prior (z0 = 1: uniform prior)
        # Cluster and local swap moves are symmetric — no correction needed
        accept = true
        if move_type == :insert
            if p_bias > 0.0
                # Composite forward density for the ACTUAL chosen site
                q_fwd = p_insert * ((insert_site in biased_set ?
                                         p_bias / length(biased_set) : 0.0) +
                                    (1.0 - p_bias) / (n_sites - n))
                ratio = z0 * p_delete / ((n + 1) * q_fwd)
            else
                # Legacy expression kept literally: bit-identical arithmetic
                ratio = z0 * (p_delete / p_insert) * (n_sites - n) / (n + 1)
            end
            if ratio < 1.0 && rand() >= ratio
                accept = false
            end
        elseif move_type == :delete
            if p_bias > 0.0
                # Reverse composite density on the POST-deletion configuration
                rev_set = lattice_biased_sites(config;
                                               predicate=bias_predicate,
                                               shells=bias_shells)
                q_rev = p_insert * ((deleted_site in rev_set ?
                                         p_bias / length(rev_set) : 0.0) +
                                    (1.0 - p_bias) / (n_sites - n + 1))
                if q_rev == 0.0
                    # Only reachable at p_bias == 1 with the vacated site
                    # outside S(x'): reject with no MH rand draw
                    accept = false
                else
                    # n is the pre-delete particle count (docstring convention)
                    ratio = n * q_rev / (z0 * p_delete)
                    if ratio < 1.0 && rand() >= ratio
                        accept = false
                    end
                end
            else
                # Legacy expression kept literally: bit-identical arithmetic
                ratio = (p_insert / p_delete) * n / (z0 * (n_sites - n + 1))
                if ratio < 1.0 && rand() >= ratio
                    accept = false
                end
            end
        end

        if accept
            # The live configuration already carries the proposal; keep it
            # (the walker's configuration field is never rebound to a copy)
            lattice.energy = proposed_energy
            n_accept += 1
            accept_this_walker = true
            if move_type == :cluster
                cluster_accepted_count += 1
            elseif move_type == :move
                swap_accepted += 1
            elseif move_type == :insert
                if insert_from_biased
                    insert_biased_accepted += 1
                else
                    insert_uniform_accepted += 1
                end
            else # :delete
                delete_accepted += 1
            end
        else
            _gc_revert_move!(config, move_type, hop_from, hop_to,
                             cluster_pairs, insert_site, deleted_site)
        end
    end

    return accept_this_walker, n_accept / max(n_steps, 1), lattice,
           cluster_accepted_count, cluster_total_count,
           (swap_attempted=swap_attempted, swap_accepted=swap_accepted,
            cluster_attempted=cluster_total_count, cluster_accepted=cluster_accepted_count,
            insert_uniform_attempted=insert_uniform_attempted,
            insert_uniform_accepted=insert_uniform_accepted,
            insert_biased_attempted=insert_biased_attempted,
            insert_biased_accepted=insert_biased_accepted,
            delete_attempted=delete_attempted, delete_accepted=delete_accepted)
end

"""
    _gc_revert_move!(config, move_type, hop_from, hop_to, cluster_pairs,
                     insert_site, deleted_site)

Internal: revert a rejected in-place grand-canonical proposal. Swap and
cluster proposals revert by involution replay (`_lattice_walk_apply!`,
`_apply_cluster_pairs!`); insertion and deletion revert by flipping the
recorded site back. No random draws.
"""
@inline function _gc_revert_move!(config, move_type::Symbol,
                                  hop_from::Int, hop_to::Int,
                                  cluster_pairs::Vector{Tuple{Int,Int}},
                                  insert_site::Int, deleted_site::Int)
    if move_type == :move
        _lattice_walk_apply!(config, hop_from, hop_to)
    elseif move_type == :cluster
        _apply_cluster_pairs!(config, cluster_pairs)
    elseif move_type == :insert
        config.components[1][insert_site] = false
    else # :delete
        config.components[1][deleted_site] = true
    end
    return config
end

"""
    gc_insert_acceptance_ratio(z0V::Float64, n::Int, p_insert::Float64, p_delete::Float64)
    gc_delete_acceptance_ratio(z0V::Float64, n::Int, p_insert::Float64, p_delete::Float64)

Metropolis acceptance ratios for the continuous-space grand-canonical kernel, with `n` the
PRE-move particle count in both directions. `z0V` is the dimensionless product of the
reference activity (dimension inverse volume) and the cell volume. The insertion ratio is
the proposal-density form z0 * p_delete / ((n + 1) * p_insert * q(r)) with the uniform
channel's q = 1/V folded in, so the shipped arithmetic is literally
z0V * (p_delete / p_insert) / (n + 1); a future biased-insertion channel enters by
evaluating a different proposal density q in place of the folded 1/V, mirroring the
composite-density pattern of the lattice kernel above.
"""
gc_insert_acceptance_ratio(z0V::Float64, n::Int, p_insert::Float64, p_delete::Float64) =
    z0V * (p_delete / p_insert) / (n + 1)

gc_delete_acceptance_ratio(z0V::Float64, n::Int, p_insert::Float64, p_delete::Float64) =
    (p_insert / p_delete) * n / z0V

"""
    MC_grand_canonical_walk!(n_steps::Int, at::AtomWalker{1}, pot::SingleComponentPotential{Pairwise}, emax::typeof(0.0u"eV");
                             z0V::Float64, species, p_move::Float64=0.5, p_insert::Float64=0.25,
                             step_size::Float64=0.5, n_max::Int=typemax(Int))

Perform a grand-canonical Monte Carlo walk on a continuous-space `AtomWalker{1}` below a
nested-sampling energy ceiling: single-atom displacements mixed with uniform-in-cell
insertions and uniform-among-particles deletions, under the Metropolis corrections that
preserve the activity-z0-weighted reference measure. The continuous counterpart
of the lattice `MC_grand_canonical_walk!` above; the site-counting acceptance ratios of the
lattice kernel do not apply here and are replaced by the volume-and-activity forms of
`gc_insert_acceptance_ratio`/`gc_delete_acceptance_ratio` (pre-move particle count in both).

Requirements: a single unfrozen component, and an orthorhombic cell (consistent with
`pbc_dist`); both are validated on entry.

Random-number stream contract (fixed by tests): one channel draw per step; a displacement
draws one particle index and the three walk displacements; an insertion draws three position
uniforms (x, y, z) plus the Metropolis uniform only when its ratio is below one; a deletion
draws one particle index plus the conditional Metropolis uniform. Guard skips (a
displacement or deletion proposed at N = 0, an insertion proposed above `n_max`) consume
only the channel draw and are not counted as attempts. The energy ceiling is checked before
the Metropolis ratio, so ceiling rejections draw no Metropolis uniform. Energies are updated
incrementally through `single_site_energy` on the same audited path the displacement walks
use, with insertions evaluated by insert-then-revert.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{1}`: The walker to evolve; a single unfrozen component.
- `pot::SingleComponentPotential{Pairwise}`: The pairwise potential.
- `emax::typeof(0.0u"eV")`: The nested-sampling energy ceiling (strict-below acceptance).
- `z0V::Float64`: Dimensionless product of the reference activity and the cell volume.
- `species`: Chemical identity (a `Symbol` or `ChemicalSpecies`) of inserted particles;
  passed explicitly because the configuration may be empty.
- `p_move::Float64=0.5`: Probability of a displacement move.
- `p_insert::Float64=0.25`: Probability of an insertion move (deletion takes the remainder).
- `step_size::Float64=0.5`: Maximum displacement per direction, in Angstrom.
- `n_max::Int=typemax(Int)`: Upper bound on the particle count, for bounded constructions
  only. A finite cap truncates the unbounded particle-number support of the reference
  measure and biases the evidence; production callers must leave it inactive.

# Returns
- `accept_this_walker::Bool`: Whether at least one move was accepted.
- `accept_rate::Float64`: Fraction of accepted moves over `n_steps`.
- `at::AtomWalker{1}`: The updated walker.
- `move_stats::NamedTuple`: Per-move-type attempt/accept counters
  (`move_*`, `insert_*`, `delete_*`); attempts are counted at proposal, ceiling rejections
  included, guard skips excluded.
"""
function MC_grand_canonical_walk!(n_steps::Int,
                                  at::AtomWalker{1},
                                  pot::SingleComponentPotential{Pairwise},
                                  emax::typeof(0.0u"eV");
                                  z0V::Float64,
                                  species::Union{Symbol, ChemicalSpecies},
                                  p_move::Float64=0.5,
                                  p_insert::Float64=0.25,
                                  step_size::Float64=0.5,
                                  n_max::Int=typemax(Int),
                                  p_bias::Float64=0.0,
                                  bias_radius::Float64=0.0,
                                  bias_grid::Int=0)
    if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
        throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
    end
    if z0V <= 0.0
        throw(ArgumentError("z0V must be positive"))
    end
    if p_bias < 0.0 || p_bias > 1.0
        throw(ArgumentError("p_bias must lie in [0, 1]"))
    end
    if p_bias > 0.0 && (bias_radius <= 0.0 || bias_grid < 1)
        throw(ArgumentError("a biased insertion channel (p_bias > 0) requires bias_radius > 0 and bias_grid >= 1"))
    end
    if p_bias == 1.0
        @warn "p_bias = 1: any deletion whose vacated cell is not a post-deletion cavity cell auto-rejects, which can freeze dense states; a mixed channel (p_bias < 1) keeps the chain irreducible." maxlog=1
    end
    if any(at.frozen)
        throw(ArgumentError("the continuous grand-canonical kernel requires an unfrozen single-component walker"))
    end
    config = at.configuration
    cellv = cell_vectors(config)
    for i in 1:3, j in 1:3
        if i != j && !iszero(ustrip(cellv[i][j]))
            throw(ArgumentError("the continuous grand-canonical kernel assumes an orthorhombic cell (consistent with pbc_dist); found a nonzero off-diagonal cell component"))
        end
    end
    box = (cellv[1][1], cellv[2][2], cellv[3][3])

    n_accept = 0
    accept_this_walker = false
    p_delete = 1.0 - p_move - p_insert
    move_attempted = 0
    move_accepted = 0
    insert_attempted = 0
    insert_accepted = 0
    insert_biased_attempted = 0
    insert_biased_accepted = 0
    delete_attempted = 0
    delete_accepted = 0

    for _ in 1:n_steps
        r = rand()
        n = at.list_num_par[1]
        if r < p_move
            # Displacement: ceiling check only (nested-sampling walk, no Metropolis)
            n == 0 && continue          # guard skip: nothing to move
            move_attempted += 1
            i_at = rand(1:n)
            prewalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
            orig_pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
            pos = single_atom_random_walk!(orig_pos, step_size)
            pos = periodic_boundary_wrap!(pos, config)
            config.position[i_at] = pos
            postwalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
            proposed_energy = at.energy + (postwalk_energy - prewalk_energy)
            if proposed_energy >= emax
                config.position[i_at] = orig_pos
            else
                at.energy = proposed_energy
                n_accept += 1
                move_accepted += 1
                accept_this_walker = true
            end
        elseif r < p_move + p_insert
            # Insertion: uniform in the cell (optionally cavity-biased),
            # insert-evaluate-revert. The p_bias = 0 path keeps the shipped
            # expressions literally and draws behind p_bias > 0 guards only.
            (n + 1 > n_max || p_insert <= 0.0) && continue   # guard skip
            insert_attempted += 1
            biased_draw = false
            q_fwd = 1.0
            if p_bias > 0.0
                cav = continuous_cavity_cells(config, box, bias_grid, bias_radius)
                biased_draw = rand() < p_bias
                if biased_draw
                    insert_biased_attempted += 1
                    if isempty(cav)
                        # Null proposal: counted as attempted, no further draws
                        continue
                    end
                    pos = _cavity_cell_draw(cav[rand(1:length(cav))], box, bias_grid)
                else
                    pos = SVector(rand() * box[1], rand() * box[2], rand() * box[3])
                end
                q_fwd = _composite_cavity_density(pos, cav, box, bias_grid, p_bias)
            else
                pos = SVector(rand() * box[1], rand() * box[2], rand() * box[3])
            end
            insert_particle!(at, pos, species)
            e_site = single_site_energy(n + 1, config, pot, at.list_num_par)
            proposed_energy = at.energy + e_site
            accept = true
            if proposed_energy >= emax
                accept = false
            else
                if p_bias > 0.0
                    ratio = gc_insert_acceptance_ratio(z0V, n, p_insert, p_delete) / q_fwd
                else
                    ratio = gc_insert_acceptance_ratio(z0V, n, p_insert, p_delete)
                end
                if ratio < 1.0 && rand() >= ratio
                    accept = false
                end
            end
            if accept
                at.energy = proposed_energy
                n_accept += 1
                insert_accepted += 1
                biased_draw && (insert_biased_accepted += 1)
                accept_this_walker = true
            else
                remove_particle!(at, n + 1)
            end
        else
            # Deletion: uniform among particles; nothing mutates until acceptance
            (n == 0 || p_delete <= 0.0) && continue          # guard skip
            delete_attempted += 1
            i_at = rand(1:n)
            e_site = single_site_energy(i_at, config, pot, at.list_num_par)
            proposed_energy = at.energy - e_site
            accept = true
            if proposed_energy >= emax
                accept = false
            else
                if p_bias > 0.0
                    # Reverse density on the post-deletion configuration; a zero
                    # reverse density rejects without a Metropolis draw
                    cav_post = continuous_cavity_cells(config, box, bias_grid, bias_radius; skip=i_at)
                    q_rev = _composite_cavity_density(position(config, i_at), cav_post,
                                                      box, bias_grid, p_bias)
                    if q_rev == 0.0
                        accept = false
                    else
                        ratio = gc_delete_acceptance_ratio(z0V, n, p_insert, p_delete) * q_rev
                        if ratio < 1.0 && rand() >= ratio
                            accept = false
                        end
                    end
                else
                    ratio = gc_delete_acceptance_ratio(z0V, n, p_insert, p_delete)
                    if ratio < 1.0 && rand() >= ratio
                        accept = false
                    end
                end
            end
            if accept
                remove_particle!(at, i_at)
                at.energy = proposed_energy
                n_accept += 1
                delete_accepted += 1
                accept_this_walker = true
            end
        end
    end

    return accept_this_walker, n_accept / max(n_steps, 1), at,
           (move_attempted=move_attempted, move_accepted=move_accepted,
            insert_attempted=insert_attempted, insert_accepted=insert_accepted,
            insert_biased_attempted=insert_biased_attempted,
            insert_biased_accepted=insert_biased_accepted,
            delete_attempted=delete_attempted, delete_accepted=delete_accepted)
end
"""
    MC_muVT_walk!(n_steps::Int, at::AtomWalker{1}, pot::SingleComponentPotential{Pairwise},
                  temperature::Float64;
                  zV::Float64, species, p_move::Float64=0.5, p_insert::Float64=0.25,
                  step_size::Float64=0.5, kb::Float64=8.617333262e-5)

Perform a fixed-temperature grand-canonical (μVT) Metropolis walk on a continuous-space
`AtomWalker{1}`: single-atom displacements under the Boltzmann acceptance
min(1, e^(-βΔU)), mixed with uniform-in-cell insertions and uniform-among-particles
deletions under min(1, ratio × e^(-βΔU)), where the ratio is the volume-and-activity
form of `gc_insert_acceptance_ratio`/`gc_delete_acceptance_ratio` (pre-move particle
count in both) and ΔU is the post-minus-pre energy change of the proposed move in every
channel. The fixed-temperature sibling of the nested-sampling kernel above, sharing its
entry validation, guard-skip discipline, incremental `single_site_energy` path,
insert-evaluate-revert bookkeeping, and `move_stats` schema; the athermal energy ceiling
is replaced by the Boltzmann factor, and the dimensionless activity-volume
`zV = e^(βμ) V / Λ(T)^3` is folded by the caller, mirroring the z0V convention (the
kernel itself never sees μ or Λ).

Draw discipline: one channel draw per step; the Metropolis uniform is drawn only when
the combined acceptance factor is below one (for displacements, only when ΔU > 0).
An insertion landing exactly on a particle (a NaN pair energy) rejects without
drawing; a +Inf overlap energy underflows the Boltzmann factor to a combined
factor of exactly zero and rejects through the ordinary draw.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{1}`: The walker to evolve; a single unfrozen component.
- `pot::SingleComponentPotential{Pairwise}`: The pairwise potential.
- `temperature::Float64`: The temperature in Kelvin; must be positive.
- `zV::Float64`: Dimensionless product of the target activity and the cell volume.
- `species`: Chemical identity (a `Symbol` or `ChemicalSpecies`) of inserted particles.
- `p_move::Float64=0.5`: Probability of a displacement move.
- `p_insert::Float64=0.25`: Probability of an insertion move (deletion takes the remainder).
- `step_size::Float64=0.5`: Maximum displacement per direction, in Angstrom.
- `kb::Float64=8.617333262e-5`: The Boltzmann constant in eV/K.

# Returns
- `at::AtomWalker{1}`: The updated walker.
- `accept_rate::Float64`: Fraction of accepted moves over `n_steps`.
- `move_stats::NamedTuple`: Per-move-type attempt/accept counters
  (`move_*`, `insert_*`, `delete_*`); attempts are counted at proposal, guard skips
  excluded.
"""
function MC_muVT_walk!(n_steps::Int,
                       at::AtomWalker{1},
                       pot::SingleComponentPotential{Pairwise},
                       temperature::Float64;
                       zV::Float64,
                       species::Union{Symbol, ChemicalSpecies},
                       p_move::Float64=0.5,
                       p_insert::Float64=0.25,
                       step_size::Float64=0.5,
                       kb::Float64=8.617333262e-5)
    if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
        throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
    end
    if zV <= 0.0
        throw(ArgumentError("zV must be positive"))
    end
    if temperature <= 0.0
        throw(ArgumentError("temperature must be positive"))
    end
    if any(at.frozen)
        throw(ArgumentError("the muVT kernel requires an unfrozen single-component walker"))
    end
    config = at.configuration
    cellv = cell_vectors(config)
    for i in 1:3, j in 1:3
        if i != j && !iszero(ustrip(cellv[i][j]))
            throw(ArgumentError("the muVT kernel assumes an orthorhombic cell (consistent with pbc_dist); found a nonzero off-diagonal cell component"))
        end
    end
    box = (cellv[1][1], cellv[2][2], cellv[3][3])

    β = 1.0 / (kb * temperature)
    n_accept = 0
    p_delete = 1.0 - p_move - p_insert
    move_attempted = 0
    move_accepted = 0
    insert_attempted = 0
    insert_accepted = 0
    delete_attempted = 0
    delete_accepted = 0

    for _ in 1:n_steps
        r = rand()
        n = at.list_num_par[1]
        if r < p_move
            # Displacement: Boltzmann acceptance; uniform drawn only for uphill moves
            n == 0 && continue          # guard skip: nothing to move
            move_attempted += 1
            i_at = rand(1:n)
            prewalk_energy = single_site_energy(i_at, config, pot, at.list_num_par)
            orig_pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
            pos = single_atom_random_walk!(orig_pos, step_size)
            pos = periodic_boundary_wrap!(pos, config)
            config.position[i_at] = pos
            dU = ustrip(u"eV", single_site_energy(i_at, config, pot, at.list_num_par) - prewalk_energy)
            if dU <= 0.0 || (isfinite(dU) && rand() < exp(-β * dU))
                at.energy = at.energy + dU * u"eV"
                n_accept += 1
                move_accepted += 1
            else
                config.position[i_at] = orig_pos
            end
        elseif r < p_move + p_insert
            # Insertion: uniform in the cell, insert-evaluate-revert
            p_insert <= 0.0 && continue   # guard skip
            insert_attempted += 1
            pos = SVector(rand() * box[1], rand() * box[2], rand() * box[3])
            insert_particle!(at, pos, species)
            e_site = ustrip(u"eV", single_site_energy(n + 1, config, pot, at.list_num_par))
            combined = gc_insert_acceptance_ratio(zV, n, p_insert, p_delete) * exp(-β * e_site)
            if combined >= 1.0 || (isfinite(combined) && rand() < combined)
                at.energy = at.energy + e_site * u"eV"
                n_accept += 1
                insert_accepted += 1
            else
                remove_particle!(at, n + 1)
            end
        else
            # Deletion: uniform among particles; nothing mutates until acceptance
            (n == 0 || p_delete <= 0.0) && continue          # guard skip
            delete_attempted += 1
            i_at = rand(1:n)
            e_site = ustrip(u"eV", single_site_energy(i_at, config, pot, at.list_num_par))
            # ΔU = -e_site, so the Boltzmann factor is e^(+β e_site)
            combined = gc_delete_acceptance_ratio(zV, n, p_insert, p_delete) * exp(β * e_site)
            if combined >= 1.0 || (isfinite(combined) && rand() < combined)
                remove_particle!(at, i_at)
                at.energy = at.energy - e_site * u"eV"
                n_accept += 1
                delete_accepted += 1
            end
        end
    end

    return at, n_accept / max(n_steps, 1),
           (move_attempted=move_attempted, move_accepted=move_accepted,
            insert_attempted=insert_attempted, insert_accepted=insert_accepted,
            delete_attempted=delete_attempted, delete_accepted=delete_accepted)
end

# Householder reflection of the 3N-space velocity off the constraint boundary:
# v' = v - 2 (v·ĝ) ĝ with ĝ the normalized energy gradient. Norm-preserving and
# an involution; both properties are pinned in the tests.
function _galilean_reflect(v::Vector{SVector{3,Float64}}, grad)
    gn2 = 0.0
    for gi in grad
        gn2 += ustrip(u"eV/Å", gi[1])^2 + ustrip(u"eV/Å", gi[2])^2 + ustrip(u"eV/Å", gi[3])^2
    end
    gn2 <= 0.0 && return v, false
    gn = sqrt(gn2)
    vg = 0.0
    for i in eachindex(v)
        vg += v[i][1] * ustrip(u"eV/Å", grad[i][1]) +
              v[i][2] * ustrip(u"eV/Å", grad[i][2]) +
              v[i][3] * ustrip(u"eV/Å", grad[i][3])
    end
    vg /= gn
    return [v[i] - (2 * vg / gn) * SVector(ustrip(u"eV/Å", grad[i][1]),
                                           ustrip(u"eV/Å", grad[i][2]),
                                           ustrip(u"eV/Å", grad[i][3])) for i in eachindex(v)], true
end

# One Galilean trajectory of n_refresh segments: straight-line motion, specular
# reflection at the first rejected trial with continuation THROUGH it, velocity
# reversal when the reflected trial also fails (continuing from the rejected
# trial, never retrying from the segment start, is what makes the segment
# retrace exactly under velocity reversal). Mutates the walker; returns the
# final velocity and the number of segments that ended inside the constraint.
function _galilean_trajectory!(at::AtomWalker{1},
                               pot::SingleComponentPotential{Pairwise},
                               emax::typeof(0.0u"eV"),
                               v::Vector{SVector{3,Float64}},
                               n_refresh::Int,
                               step_size::Float64)
    config = at.configuration
    n = at.list_num_par[1]
    n_inside = 0
    current_E = at.energy
    for _ in 1:n_refresh
        saved = copy(config.position)
        saved_E = current_E
        for i in 1:n
            pos = config.position[i] + (step_size * v[i]) * u"Å"
            config.position[i] = periodic_boundary_wrap!(pos, config)
        end
        E1 = interacting_energy(config, pot, at.list_num_par, at.frozen) + at.energy_frozen_part
        if E1 < emax
            current_E = E1
            n_inside += 1
            continue
        end
        grad = interacting_gradient(config, pot, at.list_num_par, at.frozen)
        v_r, ok = _galilean_reflect(v, grad)
        if ok
            for i in 1:n
                pos = config.position[i] + (step_size * v_r[i]) * u"Å"
                config.position[i] = periodic_boundary_wrap!(pos, config)
            end
            E2 = interacting_energy(config, pot, at.list_num_par, at.frozen) + at.energy_frozen_part
            if E2 < emax
                current_E = E2
                n_inside += 1
                v = v_r
                continue
            end
        end
        copyto!(config.position, saved)
        current_E = saved_E
        v = [-vi for vi in v]
    end
    at.energy = current_E
    return v, n_inside
end

"""
    MC_galilean_walk!(n_steps::Int, at::AtomWalker{1}, pot::SingleComponentPotential{Pairwise},
                      emax::typeof(0.0u"eV"); step_size::Float64, n_refresh::Int=8)

Perform a Galilean (reflective) Monte Carlo walk under a nested-sampling energy
ceiling: `n_steps` trajectories, each drawing a fresh velocity uniformly on the
3N-sphere and running `n_refresh` straight-line segments of 3N-space length
`step_size` (per-atom displacement scales as `step_size`/sqrt(3N)), reflecting
specularly off the constraint boundary through the energy gradient at the first
rejected trial and continuing through it, with velocity reversal when the
reflected trial also fails. All free particles move collectively; segment
energies are full `interacting_energy` recomputes (a collective move has no
incremental path), so the walker's energy never accumulates drift. The gradient
enters only through the reflection direction, so an inaccurate force degrades
mixing but cannot bias the stationary measure.

Requirements: a single unfrozen component and a fully periodic orthorhombic cell
(the position wrap mirrors coordinates at non-periodic walls without reversing
the velocity, which would break reversibility); an empty walker returns unchanged.

# Returns
- `accept_this_walker::Bool`: Whether any segment ended inside the constraint.
- `accept_rate::Float64`: Fraction of segments ending inside the constraint.
- `at::AtomWalker{1}`: The updated walker.
"""
function MC_galilean_walk!(n_steps::Int,
                           at::AtomWalker{1},
                           pot::SingleComponentPotential{Pairwise},
                           emax::typeof(0.0u"eV");
                           step_size::Float64,
                           n_refresh::Int=8)
    if step_size <= 0.0 || n_refresh <= 0
        throw(ArgumentError("step_size and n_refresh must be positive"))
    end
    if any(at.frozen)
        throw(ArgumentError("the Galilean walk kernel requires an unfrozen single-component walker"))
    end
    if !all(periodicity(at.configuration))
        throw(ArgumentError("the Galilean walk kernel requires fully periodic boundary conditions: the position wrap mirrors coordinates at non-periodic walls without reversing the velocity, which breaks the reflective walk's reversibility"))
    end
    cellv = cell_vectors(at.configuration)
    for i in 1:3, j in 1:3
        if i != j && !iszero(ustrip(cellv[i][j]))
            throw(ArgumentError("the Galilean walk kernel assumes an orthorhombic cell (consistent with pbc_dist); found a nonzero off-diagonal cell component"))
        end
    end
    n = at.list_num_par[1]
    n == 0 && return false, 0.0, at

    total_inside = 0
    for _ in 1:n_steps
        raw = [SVector(randn(), randn(), randn()) for _ in 1:n]
        nrm = sqrt(sum(vi -> vi[1]^2 + vi[2]^2 + vi[3]^2, raw))
        v = [vi / nrm for vi in raw]
        _, n_inside = _galilean_trajectory!(at, pot, emax, v, n_refresh, step_size)
        total_inside += n_inside
    end
    return total_inside > 0, total_inside / (n_steps * n_refresh), at
end

"""
    continuous_cavity_cells(config::AbstractSystem, box, grid::Int, r_cavity::Float64;
                            skip::Int=0)

The cavity cells of a continuous configuration: the sorted linear indices (of
the `grid`³ overlay, column-major) of cells whose centers clear the
minimum-image distance test r ≥ `r_cavity` against every particle. `skip`
excludes one particle index, evaluating the post-deletion cavity set without
mutating (the deletion channel's mutate-nothing-until-accepted convention).
The indicator is a deterministic function of the configuration, evaluated
identically in the biased draw, the forward proposal density, and the reverse
density; a stochastic cavity-volume estimate has no Metropolis-Hastings
exactness result and is deliberately not offered. Marks cells by particle,
O(grid³ + N (r_cavity/h)³) with h the cell edge.
"""
function continuous_cavity_cells(config::AbstractSystem, box, grid::Int, r_cavity::Float64;
                                 skip::Int=0)
    n = length(config)
    occupied = falses(grid, grid, grid)
    L = (ustrip(u"Å", box[1]), ustrip(u"Å", box[2]), ustrip(u"Å", box[3]))
    h = (L[1] / grid, L[2] / grid, L[3] / grid)
    pbc = periodicity(config)
    li = LinearIndices((grid, grid, grid))
    for p_idx in 1:n
        p_idx == skip && continue
        pos = (ustrip(u"Å", position(config, p_idx)[1]),
               ustrip(u"Å", position(config, p_idx)[2]),
               ustrip(u"Å", position(config, p_idx)[3]))
        base = (floor(Int, pos[1] / h[1]), floor(Int, pos[2] / h[2]), floor(Int, pos[3] / h[3]))
        reach = (ceil(Int, r_cavity / h[1]) + 1, ceil(Int, r_cavity / h[2]) + 1,
                 ceil(Int, r_cavity / h[3]) + 1)
        for dx in -reach[1]:reach[1], dy in -reach[2]:reach[2], dz in -reach[3]:reach[3]
            ix, iy, iz = base[1] + dx, base[2] + dy, base[3] + dz
            if pbc[1]
                ix = mod(ix, grid)
            elseif ix < 0 || ix >= grid
                continue
            end
            if pbc[2]
                iy = mod(iy, grid)
            elseif iy < 0 || iy >= grid
                continue
            end
            if pbc[3]
                iz = mod(iz, grid)
            elseif iz < 0 || iz >= grid
                continue
            end
            occupied[ix + 1, iy + 1, iz + 1] && continue
            cx, cy, cz = (ix + 0.5) * h[1], (iy + 0.5) * h[2], (iz + 0.5) * h[3]
            δx = cx - pos[1]
            pbc[1] && (δx -= L[1] * round(δx / L[1]))
            δy = cy - pos[2]
            pbc[2] && (δy -= L[2] * round(δy / L[2]))
            δz = cz - pos[3]
            pbc[3] && (δz -= L[3] * round(δz / L[3]))
            if δx^2 + δy^2 + δz^2 < r_cavity^2
                occupied[ix + 1, iy + 1, iz + 1] = true
            end
        end
    end
    return findall(!, vec(occupied))
end

# Column-major linear cell index of a position under the grid overlay, matching
# `continuous_cavity_cells` (cells located by floor on wrapped coordinates)
function _cavity_cell_index(pos, box, grid::Int)
    ix = min(floor(Int, ustrip(u"Å", pos[1]) / (ustrip(u"Å", box[1]) / grid)), grid - 1)
    iy = min(floor(Int, ustrip(u"Å", pos[2]) / (ustrip(u"Å", box[2]) / grid)), grid - 1)
    iz = min(floor(Int, ustrip(u"Å", pos[3]) / (ustrip(u"Å", box[3]) / grid)), grid - 1)
    return LinearIndices((grid, grid, grid))[max(ix, 0) + 1, max(iy, 0) + 1, max(iz, 0) + 1]
end

# Composite proposal density at a point, in units of the uniform density:
# p_bias (grid³/n_cav) 1{cell ∈ cavity} + (1 − p_bias). Membership is evaluated
# at the CELL level in both directions (a point-level test disagrees with the
# cell-level draw on boundary cells and breaks the density bookkeeping).
function _composite_cavity_density(pos, cav::Vector{Int}, box, grid::Int, p_bias::Float64)
    isempty(cav) && return 1.0 - p_bias
    in_cav = insorted(_cavity_cell_index(pos, box, grid), cav)
    return p_bias * (grid^3 / length(cav)) * (in_cav ? 1.0 : 0.0) + (1.0 - p_bias)
end

# Uniform draw within one cavity cell: the cell origin plus per-axis jitter
function _cavity_cell_draw(cell::Int, box, grid::Int)
    ci = CartesianIndices((grid, grid, grid))[cell]
    h1 = box[1] / grid
    h2 = box[2] / grid
    h3 = box[3] / grid
    return SVector((ci[1] - 1 + rand()) * h1, (ci[2] - 1 + rand()) * h2,
                   (ci[3] - 1 + rand()) * h3)
end
