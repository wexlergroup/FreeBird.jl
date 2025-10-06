"""
    MetropolisMCParameters <: SamplingParameters

Parameters for the Metropolis Monte Carlo algorithm.

# Fields
- `temperature::Float64`: The temperature of the system.
- `equilibrium_steps::Int64`: The number of steps to equilibrate the system.
- `sampling_steps::Int64`: The number of steps to sample the system.
- `step_size::Float64`: The step size for the random walk (for atomistic systems).
- `step_size_lo::Float64`: The lower bound of the step size.
- `step_size_up::Float64`: The upper bound of the step size.
- `accept_range::Tuple{Float64, Float64}`: The range of acceptance rates for adjusting the step size.
e.g. (0.25, 0.75) means that the step size will decrease if the acceptance rate is below 0.25 and increase if it is above 0.75.
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct MetropolisMCParameters <: SamplingParameters
    temperatures::Vector{Float64}
    equilibrium_steps::Int64
    sampling_steps::Int64
    step_size::Float64
    step_size_lo::Float64
    step_size_up::Float64
    accept_range::Tuple{Float64, Float64}
    random_seed::Int64
    function MetropolisMCParameters(
        temperatures;
        equilibrium_steps::Int64=10_000,
        sampling_steps::Int64=10_000,
        step_size::Float64=0.01,
        step_size_lo::Float64=0.001,
        step_size_up::Float64=1.0,
        accept_range::Tuple{Float64, Float64}=(0.25, 0.75),
        random_seed::Int64=1234
    )
        new(temperatures, equilibrium_steps, sampling_steps, step_size, step_size_lo, step_size_up, accept_range, random_seed)
    end
end


"""
    nvt_monte_carlo(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        temperature::Float64,
        num_steps::Int64,
        random_seed::Int64;
        kb::Float64 = 8.617_333_262e-5  # eV K-1
    )

Perform the NVT Monte Carlo algorithm to sample the lattice configurations.

Note: The Boltzmann constant is set to 8.617333262e-5 eV K\$^{-1}\$. Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV (defined in the Hamiltonian).

# Arguments
- `lattice::AbstractLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `random_seed::Int64`: The seed for the random number generator.
- `kb::Float64`: The Boltzmann constant (default is 8.617333262e-5 eV K-1).

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
    random_seed::Int64;
    kb::Float64 = 8.617_333_262e-5  # eV K-1
)
    # Set the random seed
    Random.seed!(random_seed)
    
    energies = Vector{Float64}(undef, num_steps)
    configurations = Vector{typeof(lattice)}(undef, num_steps)
    # df = DataFrame(energy=Float64[], config=Vector{Vector{Bool}}[])
    accepted_steps = 0

    current_lattice = deepcopy(lattice)
    current_energy = interacting_energy(current_lattice, h).val
    
    for i in 1:num_steps
        
        # Propose a swap in occupation state (only if it maintains constant N)
        proposed_lattice = deepcopy(current_lattice)

        # generate_random_new_lattice_sample!(proposed_lattice)
        lattice_random_walk!(proposed_lattice)

        # Calculate the proposed energy
        proposed_energy = interacting_energy(proposed_lattice, h).val

        # Metropolis-Hastings acceptance criterion
        # kb = 8.617_333_262e-5  # eV K-1

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
    nvt_monte_carlo(
        walker::AtomWalker,
        pot::AbstractPotential,
        temperature::Float64,
        num_steps::Int64,
        step_size::Float64,
        random_seed::Int64;
        kb::Float64 = 8.617_333_262e-5  # eV K-1
    )

Perform the NVT Monte Carlo algorithm to sample the atom walker configurations.
Note: The Boltzmann constant is set to 8.617333262e-5 eV K\$^{-1}\$. Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV (defined in the Hamiltonian).

# Arguments
- `walker::AtomWalker`: The initial atom walker configuration.
- `pot::AbstractPotential`: The potential energy function for the atoms.
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `step_size::Float64`: The step size for the random walk.
- `random_seed::Int64`: The seed for the random number generator.
- `kb::Float64`: The Boltzmann constant (default is 8.617333262e-5 eV K-1).

# Returns
- `energies::Vector{typeof(walker.energy)}`: The energies of the system at each step.
- `configurations::Vector{typeof(walker)}`: The configurations of the system at each step.
"""
function nvt_monte_carlo(
    walker::AtomWalker,
    pot::AbstractPotential,
    temperature::Float64,
    num_steps::Int64,
    step_size::Float64,
    random_seed::Int64;
    kb::Float64 = 8.617_333_262e-5  # eV K-1
)
    # Set the random seed
    Random.seed!(random_seed)

    e_unit = unit(walker.energy)
    
    energies = Vector{typeof(walker.energy)}(undef, num_steps)
    configurations = Vector{typeof(walker)}(undef, num_steps)
    accepted_steps = 0

    current_walker = deepcopy(walker)
    current_energy = interacting_energy(current_walker.configuration, pot, current_walker.list_num_par, current_walker.frozen) + current_walker.energy_frozen_part
    
    for i in 1:num_steps
        proposed_walker = deepcopy(current_walker)
        # move the atoms by 1% of the average cell size
        proposed_walker, ΔE = single_atom_random_walk!(proposed_walker, pot, step_size)
        # Calculate the proposed energy
        proposed_energy = current_energy + ΔE
        # Metropolis-Hastings acceptance criterion
        # kb = 8.617_333_262e-5  # eV K-1
        # ΔE = proposed_energy - current_energy
        if ΔE < 0*e_unit || rand() < exp(-ΔE.val / (kb * temperature))
            current_walker.configuration = proposed_walker.configuration
            current_energy = proposed_energy
            accepted_steps += 1
        end
        energies[i] = current_energy
        configurations[i] = current_walker
    end

    return energies, configurations, accepted_steps
end

"""
    monte_carlo_sampling(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        mc_params::MetropolisMCParameters;
        kb::Float64 = 8.617333262e-5 # eV/K
    )

Perform the Metropolis Monte Carlo sampling algorithm for a range of temperatures.

Note: The Boltzmann constant is set to 8.617333262e-5 eV K\$^{-1}\$. Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV (defined in the Hamiltonian).

# Arguments
- `lattice::AbstractLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `mc_params::MetropolisMCParameters`: The parameters for the Metropolis Monte Carlo algorithm.
- `kb::Float64`: The Boltzmann constant in eV/K (default is 8.617333262e-5 eV/K).

# Returns
- `energies::Vector{Float64}`: The energies of the system at each temperature.
- `configs::Vector{typeof(lattice)}`: The configurations of the system at each temperature.
- `cvs::Vector{Float64}`: The heat capacities of the system at each temperature.
- `acceptance_rates::Vector{Float64}`: The acceptance rates of the system at each temperature.
"""
function monte_carlo_sampling(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    mc_params::MetropolisMCParameters;
    kb::Float64 = 8.617333262e-5 # eV/K
)
    energies = Vector{Float64}(undef, length(mc_params.temperatures))
    cvs = Vector{Float64}(undef, length(mc_params.temperatures))
    acceptance_rates = Vector{Float64}(undef, length(mc_params.temperatures))
    configs = Vector{typeof(lattice)}(undef, length(mc_params.temperatures))

    # kb = 8.617333262e-5 # eV/K

    for (i, temp) in enumerate(mc_params.temperatures)

        # Equilibrate the lattice
        equilibration_energies, equilibration_configurations, equilibration_accepted_steps = nvt_monte_carlo(
            lattice,
            h,
            temp,
            mc_params.equilibrium_steps,
            mc_params.random_seed
        )

        equi_var = var(equilibration_energies)
        equi_mean = mean(equilibration_energies)

        equi_rate = equilibration_accepted_steps / mc_params.equilibrium_steps

        @info "Temperature: $temp K, Equilibration energy: $equi_mean, Variance: $(round(equi_var; sigdigits=4)), Acceptance rate: $(round(equi_rate; sigdigits=4))"


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
        E_var = var(sampling_energies)
        Cv = E_var / (kb * temp^2)

        # Compute the acceptance rate
        acceptance_rate = sampling_accepted_steps / mc_params.sampling_steps

        # Append the results to the DataFrame
        energies[i] = E
        cvs[i] = Cv
        acceptance_rates[i] = acceptance_rate

        configs[i] = sampling_configurations[end]
        lattice = sampling_configurations[end]

        @info "Temperature: $temp K, Energy: $E, Variance: $(round(E_var; sigdigits=4)), Cv: $(round(Cv; sigdigits=4)), Acceptance rate: $(round(acceptance_rate; sigdigits=4))"
    end

    return energies, configs, cvs, acceptance_rates
end

"""
    monte_carlo_sampling(
        at::AtomWalker,
        pot::AbstractPotential,
        mc_params::MetropolisMCParameters;
        kb::Float64 = 8.617333262e-5 # eV/K
    )

Perform the Metropolis Monte Carlo sampling algorithm for a range of temperatures.

Note: The Boltzmann constant is set to 8.617333262e-5 eV K\$^{-1}\$. Thus, the units of the temperature
should be in Kelvin, and the units of the energy should be in eV.

# Arguments
- `at::AtomWalker`: The initial atom walker configuration.
- `pot::AbstractPotential`: The potential energy function.
- `mc_params::MetropolisMCParameters`: The parameters for the Metropolis Monte Carlo algorithm.

# Returns
- `energies::Vector{Float64}`: The energies of the system at each temperature.
- `configs::Vector{typeof(at)}`: The configurations of the system at each temperature.
- `cvs::Vector{Float64}`: The heat capacities of the system at each temperature.
- `acceptance_rates::Vector{Float64}`: The acceptance rates of the system at each temperature.
"""
function monte_carlo_sampling(
    at::AtomWalker,
    pot::AbstractPotential,
    mc_params::MetropolisMCParameters;
    kb::Float64 = 8.617333262e-5 # eV/K
)
    energies = Vector{Float64}(undef, length(mc_params.temperatures))
    cvs = Vector{Float64}(undef, length(mc_params.temperatures))
    acceptance_rates = Vector{Float64}(undef, length(mc_params.temperatures))
    configs = Vector{typeof(at)}(undef, length(mc_params.temperatures))

    # kb = 8.617333262e-5 # eV/K

    for (i, temp) in enumerate(mc_params.temperatures)

        # Equilibrate the lattice
        equilibration_energies, equilibration_configurations, equilibration_accepted_steps = nvt_monte_carlo(
            at,
            pot,
            temp,
            mc_params.equilibrium_steps,
            mc_params.step_size,
            mc_params.random_seed
        )

        equi_var = var(equilibration_energies)
        equi_mean = mean(equilibration_energies)

        equi_rate = equilibration_accepted_steps / mc_params.equilibrium_steps

        adjust_step_size(mc_params, equi_rate; range=mc_params.accept_range)

        @info "Temperature: $temp K, Equilibration energy: $equi_mean, Variance: $(round(equi_var.val; sigdigits=4)), Acceptance rate: $(round(equi_rate; sigdigits=4)), Step size: $(round(mc_params.step_size; sigdigits=4))"

        # Sample the lattice
        sampling_energies, sampling_configurations, sampling_accepted_steps = nvt_monte_carlo(
            equilibration_configurations[end],
            pot,
            temp,
            mc_params.sampling_steps,
            mc_params.step_size,
            mc_params.random_seed
        )

        # Compute the heat capacity
        E = mean(sampling_energies)
        E_var = var(sampling_energies)
        Cv = E_var / (kb * temp^2)

        # Compute the acceptance rate
        acceptance_rate = sampling_accepted_steps / mc_params.sampling_steps
        adjust_step_size(mc_params, acceptance_rate; range=mc_params.accept_range)

        # Append the results to the DataFrame
        energies[i] = E.val
        cvs[i] = Cv.val
        acceptance_rates[i] = acceptance_rate

        configs[i] = sampling_configurations[end]
        at = sampling_configurations[end]

        @info "Temperature: $temp K, Energy: $E, Variance: $(round(E_var.val; sigdigits=4)), Cv: $(round(Cv.val; sigdigits=4)), Acceptance rate: $(round(acceptance_rate; sigdigits=4)), Step size: $(round(mc_params.step_size; sigdigits=4))"
    end

    ls = LJAtomWalkers(configs, pot)

    return energies, ls, cvs, acceptance_rates
end