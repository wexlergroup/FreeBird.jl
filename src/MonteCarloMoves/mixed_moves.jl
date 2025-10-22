"""
    MC_mixed_moves!(n_steps::Int, at::AtomWalker{C}, pot::AbstractPotential, step_size::Float64, emax::typeof(0.0u"eV"), freq::Vector{Int}) where C

Perform a Monte Carlo mixed move consisting of random walks and atom swaps in the `AtomWalker`.
Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The `AtomWalker` object.
- `pot::AbstractPotential`: The potential object for energy calculations.
- `step_size::Float64`: The step size for the random walk moves.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.
- `freq::Vector{Int}`: A vector containing the frequencies of random walk and swap moves. The first element is the frequency of random walk moves, and the second element is the frequency of swap moves.
Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the mixed moves.
- `at::AtomWalker{C}`: The updated `AtomWalker`.
"""
function MC_mixed_moves!(
                    n_steps::Int, 
                    at::AtomWalker{C}, 
                    pot::AbstractPotential, 
                    step_size::Float64, 
                    emax::typeof(0.0u"eV"),
                    freq::Vector{Int},
                    ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        swap_prob = freq[2] / sum(freq)
        if rand() > swap_prob
            #println("Performing a swap move")
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
        elseif rand() <= swap_prob
            #println("Performing a random walk move")
            config = at.configuration
            free_comp = free_component_index(at)
            comp1, comp2 = sample(free_comp, 2, replace=false)
            ind1 = rand(comp1)
            ind2 = rand(comp2)
            two_atoms_swap!(at, ind1, ind2)
            energy = interacting_energy(config, pot, at.list_num_par, at.frozen) + at.energy_frozen_part
            if energy >= emax
                # reject the move, revert to original position
                two_atoms_swap!(at, ind1, ind2) # swap again to go back
            else
                at.energy = energy
                # accept the move
                n_accept += 1
                accept_this_walker = true
            end
        end
    end
    return accept_this_walker, n_accept/n_steps, at
end