"""
    MetropolisMCParameters <: SamplingParameters

Parameters for the Metropolis Monte Carlo algorithm.

# Fields
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct MetropolisMCParameters <: SamplingParameters
    temperatures::Vector{Float64}
    equilibrium_steps::Int64
    sampling_steps::Int64
    random_seed::Int64
    function MetropolisMCParameters(
        temperatures;
        equilibrium_steps::Int64=10_000,
        sampling_steps::Int64=10_000,
        random_seed::Int64=1234
    )
        new(temperatures, equilibrium_steps, sampling_steps, random_seed)
    end
end


"""
    nvt_monte_carlo(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        temperature::Float64,
        num_steps::Int64,
        random_seed::Int64
    )

Perform the NVT Monte Carlo algorithm to sample the lattice configurations.

Note: The Boltzmann constant is set to 8.617333262e-5 eV K\(^{-1}\). Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV (defined in the Hamiltonian).

# Arguments
- `lattice::AbstractLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `energies::Vector{Float64}`: The energies of the system at each step.
- `configurations::Vector{typeof(lattice)}`: The configurations of the system at each step.
- `accepted_steps::Int64`: The number of accepted steps.
"""
function nvt_monte_carlo(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    temperature::Float64,
    num_steps::Int64,
    random_seed::Int64
)
    # Set the random seed
    Random.seed!(random_seed)
    
    energies = Vector{Float64}(undef, num_steps)
    configurations = Vector{typeof(lattice)}(undef, num_steps)
    # df = DataFrame(energy=Float64[], config=Vector{Vector{Bool}}[])
    accepted_steps = 0

    current_lattice = deepcopy(lattice)
    current_energy = interacting_energy(current_lattice, h).val
    
    Threads.@threads for i in 1:num_steps
        
        # Propose a swap in occupation state (only if it maintains constant N)
        proposed_lattice = deepcopy(current_lattice)

        generate_random_new_lattice_sample!(proposed_lattice)

        # Calculate the proposed energy
        proposed_energy = interacting_energy(proposed_lattice, h).val

        # Metropolis-Hastings acceptance criterion
        kb = 8.617_333_262e-5  # eV K-1

        ΔE = proposed_energy - current_energy
        if ΔE < 0 || rand() < exp(-ΔE / (kb * temperature))
            current_lattice = deepcopy(proposed_lattice)
            current_energy = proposed_energy
            accepted_steps += 1
        end

        # push!(df, (current_energy, deepcopy(current_lattice.components)))
        energies[i] = current_energy
        configurations[i] = deepcopy(current_lattice)
    end

    return energies, configurations, accepted_steps
end

"""
    monte_carlo_sampling(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        mc_params::MetropolisMCParameters
    )

Perform the Metropolis Monte Carlo sampling algorithm for a range of temperatures.

Note: The Boltzmann constant is set to 8.617333262e-5 eV K\(^{-1}\). Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV (defined in the Hamiltonian).

# Arguments
- `lattice::AbstractLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `mc_params::MetropolisMCParameters`: The parameters for the Metropolis Monte Carlo algorithm.

# Returns
- `energies::Vector{Float64}`: The energies of the system at each temperature.
- `cvs::Vector{Float64}`: The heat capacities of the system at each temperature.
- `acceptance_rates::Vector{Float64}`: The acceptance rates of the system at each temperature.
"""
function monte_carlo_sampling(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    mc_params::MetropolisMCParameters
)
    energies = Vector{Float64}(undef, length(mc_params.temperatures))
    cvs = Vector{Float64}(undef, length(mc_params.temperatures))
    acceptance_rates = Vector{Float64}(undef, length(mc_params.temperatures))

    kb = 8.617333262e-5 # eV/K

    Threads.@threads for (i, temp) in collect(enumerate(mc_params.temperatures))

        # Equilibrate the lattice
        equilibration_energies, equilibration_configurations, equilibration_accepted_steps = nvt_monte_carlo(
            lattice,
            h,
            temp,
            mc_params.equilibrium_steps,
            mc_params.random_seed
        )

        # Sample the lattice
        sampling_energies, sampling_configurations, sampling_accepted_steps = nvt_monte_carlo(
            equilibration_configurations[end],
            h,
            temp,
            mc_params.sampling_steps,
            mc_params.random_seed
        )

        # Compute the heat capacity
        E = mean(sampling_energies)
        Cv = var(sampling_energies) / (kb * temp^2)

        # Compute the acceptance rate
        acceptance_rate = sampling_accepted_steps / mc_params.sampling_steps

        # Append the results to the DataFrame
        energies[i] = E
        cvs[i] = Cv
        acceptance_rates[i] = acceptance_rate

        @info "Temperature: $temp K, Energy: $E eV, Cv: $Cv, Acceptance rate: $acceptance_rate"
    end

    return energies, cvs, acceptance_rates
end