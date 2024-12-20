"""
    nvt_monte_carlo(lattice::SLattice, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, temperature::Float64, num_steps::Int64, random_seed::Int64)

Perform a Markov chain Monte Carlo simulation of a SLattice in the NVT ensemble using the Metropolis-Hastings algorithm.

# Arguments
- `lattice::SLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `temperature::Float64`: The temperature of the system.
- `num_steps::Int64`: The number of Monte Carlo steps.
- `random_seed::Int64`: The seed for the random number generator.

# Returns
- `energies::Vector{Float64}`: A vector of the energies of the system at each step.
- `configurations::Vector{SLattice}`: A vector of the configurations of the system at each step.
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
    
    energies = Float64[]
    configurations = Vector{typeof(lattice)}()
    accepted_steps = 0

    current_lattice = deepcopy(lattice)
    current_energy = interacting_energy(current_lattice, h).val
    current_energy = interacting_energy(current_lattice, h).val
    
    push!(energies, current_energy)
    push!(configurations, deepcopy(current_lattice))

    for _ in 1:num_steps
        
        # Propose a swap in occupation state (only if it maintains constant N)
        proposed_lattice = deepcopy(current_lattice)

        generate_random_new_lattice_sample!(proposed_lattice)

        # Calculate the proposed energy
        proposed_energy = interacting_energy(proposed_lattice, h).val
        proposed_energy = interacting_energy(proposed_lattice, h).val

        # Metropolis-Hastings acceptance criterion
        k_B = 8.617_333_262e-5  # eV K-1
        ΔE = proposed_energy - current_energy
        if ΔE < 0 || rand() < exp(-ΔE / (k_B * temperature))
            current_lattice = deepcopy(proposed_lattice)
            current_energy = proposed_energy
            accepted_steps += 1
        end

        push!(energies, current_energy)
        push!(configurations, deepcopy(current_lattice))
    end

    return energies, configurations, accepted_steps
end