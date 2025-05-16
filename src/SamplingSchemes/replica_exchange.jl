"""
    ReplicaExchangeParameters

A structure to hold the parameters for the replica exchange sampling scheme.

# Fields
- `temperatures::Vector{Float64}`: The temperatures of the replicas in ascending order.
- `equilibrium_steps::Int64`: The number of steps to equilibrate each replica.
- `sampling_steps::Int64`: The number of steps to sample each replica.
- `swap_interval::Int64`: The number of steps between attempted replica-exchange moves.
- `random_seed::Int64`: The seed for the random number generator.
"""
mutable struct ReplicaExchangeParameters <: SamplingParameters
    temperatures::Vector{Float64}
    equilibrium_steps::Int64
    sampling_steps::Int64
    swap_interval::Int64
    random_seed::Int64
    function ReplicaExchangeParameters(
        temperatures;
        equilibrium_steps::Int64=10_000,
        sampling_steps::Int64=10_000,
        swap_interval::Int64=100,
        random_seed::Int64=1234
    )
        new(temperatures, equilibrium_steps, sampling_steps, swap_interval, random_seed)
    end
end

"""
    replica_exchange(
        lattice::AbstractLattice,
        h::ClassicalHamiltonian,
        re_params::ReplicaExchangeParameters
    )

Perform the replica exchange sampling scheme for a lattice system.

# Arguments
- `lattice::AbstractLattice`: The initial lattice configuration.
- `h::ClassicalHamiltonian`: The Hamiltonian containing the on-site and nearest-neighbor interaction energies.
- `re_params::ReplicaExchangeParameters`: The parameters for the replica exchange sampling scheme.

# Returns
"""
function replica_exchange(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    re_params::ReplicaExchangeParameters
)
    Random.seed!(re_params.random_seed)

    nrep = length(re_params.temperatures)
    nsteps = re_params.sampling_steps
    swapint = re_params.swap_interval

    # Per‑replica state: configurations and energies
    configs = [generate_random_new_lattice_sample!(deepcopy(lattice)) for _ in 1:nrep]
    energies = [interacting_energy(configs[i], h).val for i in 1:nrep]

    # Output containers
    energy_trajs = [Float64[] for _ in 1:nrep]
    config_trajs = [Vector{typeof(lattice)}() for _ in 1:nrep]

    # Equilibration
    print(energies)
    for (i, T) in enumerate(re_params.temperatures)
        _, cfgs, _ = nvt_monte_carlo(configs[i], h, T, re_params.equilibrium_steps, re_params.random_seed)
        configs[i] = deepcopy(cfgs[end])
        energies[i] = interacting_energy(configs[i], h).val
    end
    print(energies)

    return configs
end
