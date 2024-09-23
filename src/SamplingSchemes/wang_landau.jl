"""
    wang_landau(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, 
                num_steps::Int64, flatness_criterion::Float64, f_initial::Float64, 
                f_min::Float64, energy_bins::Vector{Float64}, random_seed::Int64)

Perform the Wang-Landau algorithm to compute the density of states for a LatticeSystem.

# Arguments
- `lattice::LatticeSystem`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
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
    h::ClassicalHamiltonian,
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
    current_energy = interacting_energy(current_lattice, h).val
    
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
            proposed_energy = interacting_energy(proposed_lattice, h).val
            
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