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

"""
    MC_random_walk!(n_steps::Int, at::AtomWalker, lj::LJParameters, step_size::Float64, emax::typeof(0.0u"eV"))

Perform a Monte Carlo random walk on the atomic/molecular system.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The walker to perform the random walk on.
- `lj::LennardJonesParametersSets`: The Lennard-Jones potential parameters.
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
                    lj::LennardJonesParametersSets, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV")
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_index = free_par_index(at)
        i_at = rand(free_index)
        pos::SVector{3, typeof(0.0u"Å")} = position(config, i_at)
        orig_pos = deepcopy(pos)
        pos = single_atom_random_walk!(pos, step_size)
        pos = periodic_boundary_wrap!(pos, config)
        config.position[i_at] = pos
        energy = interacting_energy(config, lj, at.list_num_par, at.frozen) + at.energy_frozen_part
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
    MC_random_walk!(n_steps::Int, lattice::LatticeWalker, h::LatticeGasHamiltonian, emax::Float64; energy_perturb::Float64=0.0)

Perform a Monte Carlo random walk on the lattice system.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `lattice::LatticeWalker`: The walker to perform the random walk on.
- `h::LatticeGasHamiltonian`: The lattice gas Hamiltonian.
- `emax::Float64`: The maximum energy allowed for accepting a move.
- `energy_perturb::Float64=0.0`: The energy perturbation used to make degenerate configurations distinguishable.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `lattice::LatticeWalker`: The updated walker.

"""
function MC_random_walk!(n_steps::Int,
                         lattice::LatticeWalker,
                         h::ClassicalHamiltonian,
                         emax::Float64;
                         energy_perturb::Float64=0.0,
                         )

    n_accept = 0
    accept_this_walker = false
    emax = emax * unit(lattice.energy)

    for i_mc_step in 1:n_steps
        current_lattice = lattice.configuration
        # select a random site
        site_index = rand(eachindex(current_lattice.occupations))
        # propose a swap in occupation state (only if it maintains constant N)
        proposed_lattice = deepcopy(current_lattice)

        if proposed_lattice.occupations[site_index] == 1
            vacant_sites = findall(x -> x == 0, proposed_lattice.occupations)
            if length(vacant_sites) > 0
                swap_site = rand(vacant_sites)
                proposed_lattice.occupations[site_index] = 0
                proposed_lattice.occupations[swap_site] = 1
            else
                continue
            end
        else
            occupied_sites = findall(x -> x == 1, proposed_lattice.occupations)
            if length(occupied_sites) > 0
                swap_site = rand(occupied_sites)
                proposed_lattice.occupations[site_index] = 1
                proposed_lattice.occupations[swap_site] = 0
            else
                continue
            end
        end
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

function MC_new_sample!(lattice::LatticeWalker,
                        h::LatticeGasHamiltonian,
                        emax::Float64;
                        energy_perturb::Float64=0.0,
                        )

    accept_this_walker = false
    emax = emax * unit(lattice.energy)

    current_lattice = lattice.configuration
    proposed_lattice = deepcopy(current_lattice)
    current_occupations = deepcopy(current_lattice.occupations)
    number_occupied_sites = sum(current_lattice.occupations)

    proposed_lattice.occupations = [false for i in eachindex(current_lattice.occupations)]
    for i in sample(eachindex(proposed_lattice.occupations), number_occupied_sites, replace=false)
        proposed_lattice.occupations[i] = true
    end

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