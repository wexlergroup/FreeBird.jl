include("exact_enumeration.jl")

include("nvt_monte_carlo.jl")

"""
    wang_landau(lattice::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, 
                num_steps::Int64, flatness_criterion::Float64, f_initial::Float64, 
                f_min::Float64, energy_bins::Int64, random_seed::Int64)

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
- `energy_bins::Int64`: The number of bins for the energy histogram.
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
    energy_bins::Int64,
    random_seed::Int64
)
    # Set the random seed
    Random.seed!(random_seed)
    
    # Initialize density of states and histogram
    density_of_states = ones(Float64, energy_bins)
    histogram = zeros(Int64, energy_bins)
    energies = Float64[]
    configurations = Vector{LatticeSystem}()
    
    f = f_initial
    current_lattice = deepcopy(lattice)
    current_energy = interaction_energy(current_lattice, adsorption_energy, nn_energy, nnn_energy)
    energy_min, energy_max = current_energy, current_energy

    push!(energies, current_energy)
    push!(configurations, deepcopy(current_lattice))

    # Initial energy range estimation
    for _ in 1:100
        # Propose random moves to estimate energy range
        proposed_lattice = deepcopy(current_lattice)
        site_index = rand(1:length(current_lattice.occupations))
        proposed_lattice.occupations[site_index] = 1 - proposed_lattice.occupations[site_index]
        proposed_energy = interaction_energy(proposed_lattice, adsorption_energy, nn_energy, nnn_energy)
        energy_min = min(energy_min, proposed_energy)
        energy_max = max(energy_max, proposed_energy)
    end
    
    bin_width = (energy_max - energy_min) / energy_bins
    get_bin = energy -> clamp(Int(floor((energy - energy_min) / bin_width)) + 1, 1, energy_bins)

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
            
            current_bin = get_bin(current_energy)
            proposed_bin = get_bin(proposed_energy)

            # Wang-Landau acceptance criterion
            if density_of_states[proposed_bin] == 0 || rand() < exp(density_of_states[current_bin] - density_of_states[proposed_bin])
                current_lattice = deepcopy(proposed_lattice)
                current_energy = proposed_energy
            end

            current_bin = get_bin(current_energy)
            density_of_states[current_bin] += f
            histogram[current_bin] += 1

            push!(energies, current_energy)
            push!(configurations, deepcopy(current_lattice))
        end

        # Check for flatness of the non-zero elements of the histogram
        non_zero_histogram = histogram[histogram .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > flatness_criterion * mean(non_zero_histogram)
            f /= 2
            histogram .= 0
        end

        # Print progress
        println("f = $f")
    end

    # Calculate the energy for each bin
    bin_energies = [energy_min + (i - 0.5) * bin_width for i in 1:energy_bins]
    return density_of_states, histogram, bin_energies, energies, configurations
end