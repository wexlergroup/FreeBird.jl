"""
    rejection_sampling(walker::SLattice, h::ClassicalHamiltonian)

Perform rejection sampling to generate a new configuration with energy below a specified limit.

# Arguments
- `walker::SLattice`: The walker to sample from.
- `energy_limit::Float64`: The energy limit.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.

# Returns
- `walker::SLattice`: The updated walker.
- `current_energy::Float64`: The energy of the updated walker.

"""

k_B = 8.617_333_262e-5

function rejection_sampling(
    walker::SLattice, 
    energy_limit::Float64, 
    h::ClassicalHamiltonian,
    perturbation::Float64
)
    current_energy = interacting_energy(walker, h).val + 1 / perturbation

    while current_energy > energy_limit
        # Save the current lattice state
        current_occupations = deepcopy(walker.occupations)
        number_occupied_sites = sum(walker.occupations)  # Number of occupied sites

        # Generate a new walker with a random configuration
        walker.occupations = [false for i in 1:length(walker.occupations)]
        for i in sample(1:length(walker.occupations), number_occupied_sites, replace=false)
            walker.occupations[i] = true
        end
        
        # Calculate the new energy
        current_energy = interacting_energy(walker, h).val
        current_energy += perturbation * (rand() - 0.5)
        
        # If the new energy is below the energy limit, we keep this configuration
        if current_energy < energy_limit
            break
        else
            # Otherwise, revert the change and try again
            walker.occupations = deepcopy(current_occupations)
        end
    end
    return walker, current_energy
end


function monte_carlo_displacement_step_constant_N!(lattice::SLattice, h::ClassicalHamiltonian, temperature::Float64)

    # Randomly select any site
    site_index = rand(1:length(lattice.occupations))
    
    # Propose a swap in occupation state while maintaining constant N
    proposed_lattice = deepcopy(lattice)
    
    if proposed_lattice.occupations[site_index]
        # Current site is occupied, find a vacant site to swap with
        vacant_sites = findall(x -> x == false, proposed_lattice.occupations)
        if isempty(vacant_sites)
            return interacting_energy(lattice, h).val  # No change possible
        end
        swap_site = rand(vacant_sites)
        proposed_lattice.occupations[site_index] = false
        proposed_lattice.occupations[swap_site] = true
    else
        # Current site is vacant, find an occupied site to swap with
        occupied_sites = findall(x -> x == true, proposed_lattice.occupations)
        if isempty(occupied_sites)
            return interacting_energy(lattice, h).val  # No change possible
        end
        swap_site = rand(occupied_sites)
        proposed_lattice.occupations[site_index] = true
        proposed_lattice.occupations[swap_site] = false
    end

    # Calculate energy difference
    current_energy = interacting_energy(lattice, h).val
    proposed_energy = interacting_energy(proposed_lattice, h).val
    ΔE = proposed_energy - current_energy

    # Metropolis criterion
    if ΔE < 0 || exp(-ΔE / (k_B * temperature)) > rand()
        # Accept the swap
        lattice.occupations[site_index], lattice.occupations[swap_site] = proposed_lattice.occupations[site_index], proposed_lattice.occupations[swap_site]
        return proposed_energy
    end
    return current_energy
end

# Function to attempt swapping configurations between two replicas
function attempt_swap!(replicas::Vector{L}, energies::Vector{Float64}, temperatures::Vector{Float64}, index1::Int, index2::Int) where L

    Δ = (1/(k_B * temperatures[index1]) - 1/(k_B * temperatures[index2])) * (energies[index2] - energies[index1])
    if Δ < 0 || exp(-Δ) > rand()
        replicas[index1], replicas[index2] = replicas[index2], replicas[index1]
        energies[index1], energies[index2] = energies[index2], energies[index1]
        return true
    end
    return false
end

function nvt_replica_exchange(replicas::Vector{L}, temperatures::Vector{Float64}, steps::Int, h::ClassicalHamiltonian, swap_fraction::Float64) where L
    num_replicas = length(replicas)
    energies = [interacting_energy(replica, h).val for replica in replicas]
    
    # Initialize array to track which replica is at each temperature
    temperature_indices = collect(1:num_replicas)
    
    # Array to store temperature history for each replica
    temperature_history = zeros(Int, steps, num_replicas)
    
    for step in 1:steps
        for i in 1:num_replicas
            temperature_history[step, temperature_indices[i]] = temperatures[i]
            
            if rand() < swap_fraction && i < num_replicas
                # Attempt swap with next replica with a probability of swap_fraction
                success = attempt_swap!(replicas, energies, temperatures, temperature_indices[i], temperature_indices[i+1])
                if success
                    # Swap temperature indices if the swap was successful
                    temperature_indices[i], temperature_indices[i+1] = temperature_indices[i+1], temperature_indices[i]
                end
            else
                # Perform a Monte Carlo displacement step
                energies[temperature_indices[i]] = monte_carlo_displacement_step_constant_N!(replicas[temperature_indices[i]], h, temperatures[i])
            end
        end
    end
    
    return replicas, energies, temperature_history
end