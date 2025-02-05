"""
    WangLandauParameters

A structure to hold the parameters for the Wang-Landau sampling scheme.

# Fields
- `num_steps::Int64`: The number of Monte Carlo steps.
- `flatness_criterion::Float64`: The criterion for flatness of the histogram.
- `f_initial::Float64`: The initial modification factor.
- `f_min::Float64`: The minimum modification factor.
- `energy_bins::Vector{Float64}`: The pre-supplied energy bins.
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct WangLandauParameters <: SamplingParameters
    num_steps::Int64
    flatness_criterion::Float64
    f_initial::Float64
    f_min::Float64
    energy_bins::Vector{Float64}
    random_seed::Int64
    function WangLandauParameters(;
        num_steps::Int64=100,
        flatness_criterion::Float64=0.8,
        f_initial::Float64=Float64(MathConstants.e),
        f_min::Float64=exp(1e-8),
        energy_min::Float64=0.0,
        energy_max::Float64=1.0,
        num_energy_bins::Int64=100,
        random_seed::Int64=1234,
    )
        energy_bins = collect(range(energy_min, stop=energy_max, length=num_energy_bins))
        new(num_steps, flatness_criterion, f_initial, f_min, energy_bins, random_seed)   
    end
end


# Function to determine the bin index for a given energy
function get_bin_index(energy::Float64, bins::Vector{Float64})
    for i in 1:length(bins)-1
        if energy >= bins[i] && energy < bins[i+1]
            return i
        end
    end
    return length(bins)
end


"""
    wang_landau(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        wl_params::WangLandauParameters
    )

    wang_landau(
        walker::AtomWalker,
        lj::LennardJonesParametersSets,
        wl_params::WangLandauParameters
    )

Perform the Wang-Landau sampling scheme for a lattice or an atomistic system.

# Arguments
- `lattice::AbstractLattice`/`walker::AtomWalker`: The initial lattice/atomistic configuration.
- `h::ClassicalHamiltonian`/`lj::LennardJonesParametersSets`: The Hamiltonian parameters for the lattice/atomistic system.
- `wl_params::WangLandauParameters`: The parameters for the Wang-Landau sampling scheme.

# Returns
- `df::DataFrame`/`energies::Vector{Float64}`: The energies of the system at each step.
- `wl_params::WangLandauParameters`: The parameters for the Wang-Landau sampling scheme.
- `S::Vector{Float64}`: The entropy of the system.
- `H::Vector{Int64}`: The histogram of the system.
"""
function wang_landau(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    wl_params::WangLandauParameters
)
    # Set the random seed
    Random.seed!(wl_params.random_seed)
    
    energy_bins_count = length(wl_params.energy_bins)
    
    # Set g(E) = 1 and H(E) = 0
    # g = ones(Float64, energy_bins_count)
    S = zeros(Float64, energy_bins_count)  # Entropy
    H = zeros(Int64, energy_bins_count)

    df = DataFrame(energy=Float64[], config=Vector{Vector{Bool}}[])
    
    # Choose a modification factor
    f = wl_params.f_initial
    
    current_lattice = deepcopy(lattice)
    current_energy = interacting_energy(current_lattice, h).val
    
    # push!(energies, current_energy)
    # push!(configurations, deepcopy(current_lattice))
    push!(df, (current_energy, deepcopy(current_lattice.components)))
    counter = 0

    while f > wl_params.f_min

        counter += 1

        for _ in 1:wl_params.num_steps        
            # Propose a swap in occupation state (only if it maintains constant N)
            proposed_lattice = deepcopy(current_lattice)

            generate_random_new_lattice_sample!(proposed_lattice)

            # Calculate the proposed energy
            proposed_energy = interacting_energy(proposed_lattice, h).val
            
            current_bin = get_bin_index(current_energy, wl_params.energy_bins)
            proposed_bin = get_bin_index(proposed_energy, wl_params.energy_bins)

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
            current_bin = get_bin_index(current_energy, wl_params.energy_bins)
            # g[current_bin] *= f
            S[current_bin] += log(f)
            H[current_bin] += 1

            push!(df, (current_energy, deepcopy(current_lattice.components)))
        end

        # If the histogram is flat, decrease f, e.g. f_{i + 1} = f_i^{1/2}
        non_zero_histogram = H[H .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > wl_params.flatness_criterion * mean(non_zero_histogram)
            f = sqrt(f)
            H .= 0
            # Print progress
            @info "f = $f, iterations = $counter"
            counter = 0
        end
    end

    return df, wl_params, S, H
end


function wang_landau(
    walker::AtomWalker,
    lj::LennardJonesParametersSets,
    wl_params::WangLandauParameters
    )

    # Set the random seed
    Random.seed!(wl_params.random_seed)
    
    energy_bins_count = length(wl_params.energy_bins)
    
    # Set g(E) = 1 and H(E) = 0
    # g = ones(Float64, energy_bins_count)
    S = zeros(Float64, energy_bins_count)  # Entropy
    H = zeros(Int64, energy_bins_count)

    energies = Float64[]
    configs = AtomWalker[]
    
    # Choose a modification factor
    f = wl_params.f_initial

    current_walker = deepcopy(walker)
    current_energy = interacting_energy(current_walker.configuration, lj).val

    push!(energies, current_energy)
    push!(configs, deepcopy(current_walker))
    counter = 0

    while f > wl_params.f_min

        counter += 1

        for _ in 1:wl_params.num_steps        
            # Propose a swap in occupation state (only if it maintains constant N)
            proposed_walker = deepcopy(current_walker)

            MC_random_walk!(200, proposed_walker, lj, 0.01, Inf*unit(proposed_walker.energy))

            # Calculate the proposed energy
            proposed_energy = interacting_energy(proposed_walker.configuration, lj).val
            
            current_bin = get_bin_index(current_energy, wl_params.energy_bins)
            proposed_bin = get_bin_index(proposed_energy, wl_params.energy_bins)

            # Calculate the ratio of the density of states which results if the occupation state is swapped
            # η = g[current_bin] / g[proposed_bin]
            η = exp(S[current_bin] - S[proposed_bin])

            # Generate a random number r such that 0 < r < 1
            r = rand()

            # If r < η, swap the occupation state
            if r < η
                current_walker = deepcopy(proposed_walker)
                current_energy = proposed_energy
            end

            # Set g(E) = g(E) * f and H(E) = H(E) + 1
            current_bin = get_bin_index(current_energy, wl_params.energy_bins)
            # g[current_bin] *= f
            S[current_bin] += log(f)
            H[current_bin] += 1

            push!(energies, current_energy)
            push!(configs, deepcopy(current_walker))
        end

        # If the histogram is flat, decrease f, e.g. f_{i + 1} = f_i^{1/2}
        non_zero_histogram = H[H .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > wl_params.flatness_criterion * mean(non_zero_histogram)
            f = sqrt(f)
            H .= 0
            # Print progress
            @info "f = $f, iterations = $counter"
            counter = 0
        end
    end

    return energies, configs, wl_params, S, H
end