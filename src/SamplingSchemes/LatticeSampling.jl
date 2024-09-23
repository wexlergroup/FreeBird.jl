include("exact_enumeration.jl")

include("nvt_monte_carlo.jl")

include("wang_landau.jl")

"""
    rejection_sampling(walker::LatticeSystem, energy_limit::Float64, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64)

Perform rejection sampling to generate a new configuration with energy below a specified limit.

# Arguments
- `walker::LatticeSystem`: The walker to sample from.
- `energy_limit::Float64`: The energy limit.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.

# Returns
- `walker::LatticeSystem`: The updated walker.
- `current_energy::Float64`: The energy of the updated walker.

"""

function rejection_sampling(
    walker::LatticeSystem, 
    energy_limit::Float64, 
    adsorption_energy::Float64, 
    nn_energy::Float64, 
    nnn_energy::Float64,
    perturbation::Float64
)
    current_energy = interaction_energy(walker, adsorption_energy, nn_energy, nnn_energy) + 1 / perturbation

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
        current_energy = interaction_energy(walker, adsorption_energy, nn_energy, nnn_energy)
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


function monte_carlo_displacement_step_constant_N!(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, temperature::Float64)
    # Randomly select any site
    site_index = rand(1:length(lattice.occupations))
    
    # Propose a swap in occupation state while maintaining constant N
    proposed_lattice = deepcopy(lattice)
    
    if proposed_lattice.occupations[site_index]
        # Current site is occupied, find a vacant site to swap with
        vacant_sites = findall(x -> x == false, proposed_lattice.occupations)
        if isempty(vacant_sites)
            return interaction_energy(lattice, adsorption_energy, nn_energy, nnn_energy)  # No change possible
        end
        swap_site = rand(vacant_sites)
        proposed_lattice.occupations[site_index] = false
        proposed_lattice.occupations[swap_site] = true
    else
        # Current site is vacant, find an occupied site to swap with
        occupied_sites = findall(x -> x == true, proposed_lattice.occupations)
        if isempty(occupied_sites)
            return interaction_energy(lattice, adsorption_energy, nn_energy, nnn_energy)  # No change possible
        end
        swap_site = rand(occupied_sites)
        proposed_lattice.occupations[site_index] = true
        proposed_lattice.occupations[swap_site] = false
    end

    # Calculate energy difference
    current_energy = interaction_energy(lattice, adsorption_energy, nn_energy, nnn_energy)
    proposed_energy = interaction_energy(proposed_lattice, adsorption_energy, nn_energy, nnn_energy)
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
function attempt_swap!(replicas::Vector{LatticeSystem}, energies::Vector{Float64}, temperatures::Vector{Float64}, index1::Int, index2::Int)
    Δ = (1/(k_B * temperatures[index1]) - 1/(k_B * temperatures[index2])) * (energies[index2] - energies[index1])
    if Δ < 0 || exp(-Δ) > rand()
        replicas[index1], replicas[index2] = replicas[index2], replicas[index1]
        energies[index1], energies[index2] = energies[index2], energies[index1]
        return true
    end
    return false
end

function nvt_replica_exchange(replicas::Vector{LatticeSystem}, temperatures::Vector{Float64}, steps::Int, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, swap_fraction::Float64)
    num_replicas = length(replicas)
    energies = [interaction_energy(replica, adsorption_energy, nn_energy, nnn_energy) for replica in replicas]
    
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
                energies[temperature_indices[i]] = monte_carlo_displacement_step_constant_N!(replicas[temperature_indices[i]], adsorption_energy, nn_energy, nnn_energy, temperatures[i])
            end
        end
    end
    
    return replicas, energies, temperature_history
end