include("exact_enumeration.jl")

include("nvt_monte_carlo.jl")

"""
    wang_landau(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, 
                num_steps::Int64, flatness_criterion::Float64, f_initial::Float64, 
                f_min::Float64, energy_bins::Vector{Float64}, random_seed::Int64)

Perform the Wang-Landau algorithm to compute the density of states for a LatticeSystem.

# Arguments
- `lattice::LatticeSystem`: The initial lattice configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `flatness_criterion::Float64`: The criterion for flatness of the histogram.
- `f_initial::Float64`: The initial modification factor.
- `f_min::Float64`: The minimum modification factor.
- `energy_bins::Vector{Float64}`: The pre-supplied energy bins.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `density_of_states::Vector{Float64}`: The estimated density of states.
- `histogram::Vector{Int64}`: The histogram of visited energy states.
- `energies::Vector{Float64}`: The energies of the system at each step.
- `configurations::Vector{LatticeSystem}`: The configurations of the system at each step.
"""

function wang_landau(
    lattice::LatticeSystem,
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    num_steps::Int64,
    flatness_criterion::Float64,
    f_initial::Float64,
    f_min::Float64,
    energy_bins::Vector{Float64},
    random_seed::Int64
)
    # Set the random seed
    Random.seed!(random_seed)
    
    energy_bins_count = length(energy_bins)
    
    # Set g(E) = 1 and H(E) = 0
    # g = ones(Float64, energy_bins_count)
    S = zeros(Float64, energy_bins_count)  # Entropy
    H = zeros(Int64, energy_bins_count)

    energies = Float64[]
    configurations = Vector{LatticeSystem}()
    
    # Choose a modification factor
    f = f_initial
    
    current_lattice = deepcopy(lattice)
    current_energy = interaction_energy(current_lattice, adsorption_energy, nn_energy, nnn_energy)
    
    push!(energies, current_energy)
    push!(configurations, deepcopy(current_lattice))
    
    # Function to determine the bin index for a given energy
    function get_bin_index(energy::Float64, bins::Vector{Float64})
        for i in 1:length(bins)-1
            if energy >= bins[i] && energy < bins[i+1]
                return i
            end
        end
        return length(bins)
    end

    while f > f_min
        for _ in 1:num_steps

            # Choose a site
            site_index = rand(1:length(current_lattice.occupations))
            
            # Propose a swap in occupation state (only if it maintains constant N)
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

            # Calculate the proposed energy
            proposed_energy = interaction_energy(proposed_lattice, adsorption_energy, nn_energy, nnn_energy)
            
            current_bin = get_bin_index(current_energy, energy_bins)
            proposed_bin = get_bin_index(proposed_energy, energy_bins)

            # Calculate the ratio of the density of states which results if the occupation state is swapped
            # η = g[current_bin] / g[proposed_bin]
            η = exp(S[current_bin] - S[proposed_bin])

            # Generate a random number r such that 0 < r < 1
            r = rand()

            # If r < η, swap the occupation state
            if r < η
                current_lattice = deepcopy(proposed_lattice)
                current_energy = proposed_energy
            end

            # Set g(E) = g(E) * f and H(E) = H(E) + 1
            current_bin = get_bin_index(current_energy, energy_bins)
            # g[current_bin] *= f
            S[current_bin] += log(f)
            H[current_bin] += 1

            push!(energies, current_energy)
            push!(configurations, deepcopy(current_lattice))
        end

        # If the histogram is flat, decrease f, e.g. f_{i + 1} = f_i^{1/2}
        non_zero_histogram = H[H .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > flatness_criterion * mean(non_zero_histogram)
            f = sqrt(f)
            H .= 0
        end

        # Print progress
        println("f = $f")
    end

    return S, H, energy_bins, energies, configurations
end

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

"""
    nested_sampling(walkers::Vector{LatticeSystem}, num_steps::Int, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64)

Perform the nested sampling algorithm to sample the energy landscape of a LatticeSystem.

# Arguments
- `walkers::Vector{LatticeSystem}`: An array of initialized walkers.
- `num_steps::Int`: The number of steps to perform.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.

# Returns
- `walkers::Vector{LatticeSystem}`: The updated walkers.
- `energies::Vector{Float64}`: The energies of the walkers.
- `all_energies::Vector{Float64}`: The energies of all sampled configurations.

"""

function nested_sampling(
    walkers::Vector{LatticeSystem},  # User supplies an array of initialized walkers
    num_steps::Int,
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    perturbation::Float64
)
    # Calculate initial energies of the supplied walkers
    energies = [interaction_energy(walker, adsorption_energy, nn_energy, nnn_energy) + perturbation * (rand() - 0.5) for walker in walkers]

    # Array to store all sampled energies
    all_energies = Float64[]
    
    # Main nested sampling loop
    for i in 1:num_steps
        # Find the index of the walker with the highest energy
        high_energy_index = argmax(energies)
        high_energy_walker = deepcopy(walkers[high_energy_index])

        # Set the energy limit as the current highest energy
        energy_limit = energies[high_energy_index]

        # Save the highest energy
        push!(all_energies, energy_limit)

        # Generate a new walker using rejection sampling
        new_walker, new_energy = rejection_sampling(high_energy_walker, energy_limit, adsorption_energy, nn_energy, nnn_energy, perturbation)
        
        # Replace the highest energy walker if the new energy is less than the highest energy
        if new_energy < energies[high_energy_index]
            walkers[high_energy_index] = new_walker
            energies[high_energy_index] = new_energy
        end
    end
    
    return walkers, energies, all_energies
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