"""
    two_atoms_swap!(at::AtomWalker{C}, ind1, ind2) where C

Swap the positions of two atoms in the `AtomWalker`.

# Arguments
- `at::AtomWalker{C}`: The `AtomWalker` object.
- `ind1::Int`: The index of the first atom.
- `ind2::Int`: The index of the second atom.

# Returns
- `at::AtomWalker{C}`: The updated `AtomWalker`.
"""
function two_atoms_swap!(at::AtomWalker{C}, ind1, ind2) where C
    config = at.configuration
    config.position[ind1], config.position[ind2] = config.position[ind2], config.position[ind1]
    config.species[ind1], config.species[ind2] = config.species[ind2], config.species[ind1]
    return at
end

"""
    MC_random_swap!(n_steps::Int, at::AtomWalker{C}, lj::LennardJonesParametersSets, emax::typeof(0.0u"eV"))

Perform a Monte Carlo random swap of two atoms in the `AtomWalker`. Only works when there are two or more non-frozen components.

# Arguments
- `n_steps::Int`: The number of Monte Carlo steps to perform.
- `at::AtomWalker{C}`: The `AtomWalker` object.
- `lj::LennardJonesParametersSets`: The Lennard-Jones parameters.
- `emax::typeof(0.0u"eV")`: The maximum energy allowed for accepting a move.

# Returns
- `accept_this_walker::Bool`: Whether the walker is accepted or not.
- `accept_rate::Float64`: The acceptance rate of the random walk.
- `at::AtomWalker{C}`: The updated `AtomWalker`.
"""
function MC_random_swap!(n_steps::Int, 
                         at::AtomWalker{C}, 
                         lj::LennardJonesParametersSets, 
                         emax::typeof(0.0u"eV")
                         ) where C
    n_accept = 0
    accept_this_walker = false
    for i_mc_step in 1:n_steps
        config = at.configuration
        free_comp = free_component_index(at)
        comp1, comp2 = sample(free_comp, 2, replace=false)
        ind1 = rand(comp1)
        ind2 = rand(comp2)
        two_atoms_swap!(at, ind1, ind2)
        energy = interacting_energy(config, lj, at.list_num_par, at.frozen) + at.energy_frozen_part
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
    return accept_this_walker, n_accept/n_steps, at
end