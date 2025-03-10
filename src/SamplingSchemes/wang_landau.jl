"""
    WangLandauParameters

A structure to hold the parameters for the Wang-Landau sampling scheme.

# Fields
- `num_steps::Int64`: The number of Monte Carlo steps.
- `flatness_criterion::Float64`: The criterion for flatness of the histogram.
- `f_initial::Float64`: The initial modification factor.
- `f_min::Float64`: The minimum modification factor.
- `energy_bins::Vector{Float64}`: The pre-supplied energy bins.
- `max_iter::Int64`: The maximum number of iterations in each flatness check.
- `step_size::Float64`: The step size for the random walk (for atomistic systems).
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct WangLandauParameters <: SamplingParameters
    num_steps::Int64
    flatness_criterion::Float64
    f_initial::Float64
    f_min::Float64
    energy_bins::Vector{Float64}
    max_iter::Int64
    step_size::Float64
    random_seed::Int64
end

"""
    WangLandauParameters(; 
        num_steps::Int64=100,
        flatness_criterion::Float64=0.8,
        f_initial::Float64=Float64(MathConstants.e),
        f_min::Float64=exp(1e-8),
        energy_min::Float64=0.0,
        energy_max::Float64=1.0,
        num_energy_bins::Int64=100,
        max_iter::Int64=1000,
        step_size::Float64=0.01,
        random_seed::Int64=1234
    )

Create a `WangLandauParameters` object with the specified parameters.

# Arguments
- `num_steps::Int64`: The number of Monte Carlo steps.
- `flatness_criterion::Float64`: The criterion for flatness of the histogram.
- `f_initial::Float64`: The initial modification factor.
- `f_min::Float64`: The minimum modification factor.
- `energy_min::Float64`: The minimum energy.
- `energy_max::Float64`: The maximum energy.
- `num_energy_bins::Int64`: The number of energy bins.
- `max_iter::Int64`: The maximum number of iterations in each flatness check.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `WangLandauParameters`: The parameters for the Wang-Landau sampling scheme.
"""
function WangLandauParameters(;
            num_steps::Int64=100,
            flatness_criterion::Float64=0.8,
            f_initial::Float64=Float64(MathConstants.e),
            f_min::Float64=exp(1e-8),
            energy_min::Float64=0.0,
            energy_max::Float64=1.0,
            num_energy_bins::Int64=100,
            max_iter::Int64=1000,
            step_size::Float64=0.01,
            random_seed::Int64=1234,
            )
    energy_bins = collect(range(energy_min, stop=energy_max, length=num_energy_bins))
    WangLandauParameters(num_steps, flatness_criterion, f_initial, f_min, energy_bins, max_iter, step_size, random_seed)
end


# Function to determine the bin index for a given energy
function get_bin_index(energy::Float64, bins::Vector{Float64})
    return searchsortedfirst(bins, energy)
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
- `configs::Vector{AbstractLattice}`/`configs::Vector{AtomWalker}`: The configurations of the system at each step.
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
    configs = AbstractLattice[]
    
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
            push!(configs, deepcopy(current_lattice))
        end

        # If the histogram is flat, decrease f, e.g. f_{i + 1} = f_i^{1/2}
        non_zero_histogram = H[H .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > wl_params.flatness_criterion * mean(non_zero_histogram)
            f = sqrt(f)
            H .= 0
            # Print progress
            @info "f = $f, iterations = $counter"
            counter = 0
        elseif counter > wl_params.max_iter
            @warn "Maximum number of iterations reached!"
            break
        end
    end

    return df, configs, wl_params, S, H
end


function wang_landau(
    walker::AtomWalker,
    lj::LennardJonesParametersSets,
    wl_params::WangLandauParameters
    )

    # Set the random seed
    Random.seed!(wl_params.random_seed)

    # e_unit = unit(walker.energy)
    # cell_vec = walker.configuration.cell.cell_vectors
    # cell_volume = cell_vec[1][1] * cell_vec[2][2] * cell_vec[3][3]
    # cell_size = cbrt(cell_volume).val
    
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
    current_energy = interacting_energy(current_walker.configuration, lj, current_walker.list_num_par, current_walker.frozen) + current_walker.energy_frozen_part
    current_walker.energy = current_energy
    println("Initial energy: ", current_energy)
    push!(energies, current_energy.val)
    push!(configs, deepcopy(current_walker))
    counter = 0

    while f > wl_params.f_min

        counter += 1

        accepted = 0

        for _ in 1:wl_params.num_steps
            
            # Propose a swap in occupation state (only if it maintains constant N)
            proposed_walker = deepcopy(current_walker)

            proposed_walker, ΔE = single_atom_random_walk!(proposed_walker, lj, wl_params.step_size)

            # Calculate the proposed energy
            # proposed_energy = interacting_energy(proposed_walker.configuration, lj, proposed_walker.list_num_par, proposed_walker.frozen) + proposed_walker.energy_frozen_part
            # proposed_walker.energy = proposed_energy
            proposed_energy = current_energy + ΔE
            proposed_walker.energy = proposed_energy
            # @info "Proposed energy: $proposed_energy"

            if proposed_energy.val > wl_params.energy_bins[end] || proposed_energy.val < wl_params.energy_bins[1]
                continue
            end
            
            
            current_bin = get_bin_index(current_energy.val, wl_params.energy_bins)
            proposed_bin = get_bin_index(proposed_energy.val, wl_params.energy_bins)

            if current_bin == proposed_bin
                continue
            end

            # Calculate the ratio of the density of states which results if the occupation state is swapped
            # η = g[current_bin] / g[proposed_bin]
            η = exp(S[current_bin] - S[proposed_bin])

            # Generate a random number r such that 0 < r < 1
            r = rand()

            # If r < η, swap the occupation state
            if r < η
                current_walker = deepcopy(proposed_walker)
                current_energy = proposed_energy
                accepted += 1
                # push!(configs, current_walker)
                # @info "Accepted, energy = $current_energy"
            end

            # Set g(E) = g(E) * f and H(E) = H(E) + 1
            current_bin = get_bin_index(current_energy.val, wl_params.energy_bins)
            # g[current_bin] *= f
            S[current_bin] += log(f)
            H[current_bin] += 1

            push!(energies, current_energy.val)
        end

        @info "Iteration: $counter, f = $f, current_energy = $current_energy, acceptance rate = $(accepted/wl_params.num_steps)"
        
        # Save the last configuration
        push!(configs, deepcopy(current_walker)) 

        # If the histogram is flat, decrease f, e.g. f_{i + 1} = f_i^{1/2}
        non_zero_histogram = H[H .> 0]
        if length(non_zero_histogram) > 0 && minimum(non_zero_histogram) > wl_params.flatness_criterion * mean(non_zero_histogram)
            f = sqrt(f)
            # Print progress
            @info "Update f = $f, iterations = $counter"
            H .= 0
            counter = 0
        elseif counter > wl_params.max_iter
            @warn "Maximum number of iterations reached!"
            break
        end
    end

    return energies, configs, wl_params, S, H
end