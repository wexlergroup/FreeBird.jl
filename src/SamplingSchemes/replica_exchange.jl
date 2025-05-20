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
    random_seed::Int64  # TODO: pass an AbstractRNG instead of an Int64
    function ReplicaExchangeParameters(
        temperatures;
        equilibrium_steps::Int64=10_000,
        sampling_steps::Int64=10_000,
        swap_interval::Int64=100,
        random_seed::Int64=1234  # TODO: pass an AbstractRNG instead of an Int64
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
# If `store_trajectories` is `true`, returns `(energy_trajs, config_trajs,
# assignment_trajs, swap_acc_rate)` where `assignment_trajs` contains the
# temperature index associated with each saved configuration.
# Otherwise returns `(energies, configs, swap_acc_rate)`.
"""
function replica_exchange(
    lattice::AbstractLattice,
    h::ClassicalHamiltonian,
    re_params::ReplicaExchangeParameters;
    store_trajectories::Bool = true
)

    Random.seed!(re_params.random_seed)

    nrep = length(re_params.temperatures)
    nsteps = re_params.sampling_steps
    swapint = re_params.swap_interval
    println("Replica exchange sampling with $nrep replicas, $nsteps steps, and swap interval of $swapint.")

    # Per‑replica state: configurations and energies
    configs = [generate_random_new_lattice_sample!(deepcopy(lattice)) for _ in 1:nrep]
    energies = [interacting_energy(configs[i], h).val for i in 1:nrep]
    config_ids = collect(1:nrep)  # Track identity of each configuration
    println("Initial energies: ", energies)

    # Output containers
    energy_trajs = store_trajectories ? [Float64[] for _ in 1:nrep] : nothing
    config_trajs = store_trajectories ? [Vector{typeof(lattice)}() for _ in 1:nrep] : nothing
    assignment_trajs = store_trajectories ? [Int[] for _ in 1:nrep] : nothing
    println(
        "Output containers: ",
        size(energy_trajs), " ", size(config_trajs), " ", size(assignment_trajs)
    )

    # ──────────────────────── Equilibration ────────────────────────
    for (i, T) in enumerate(re_params.temperatures)
        seed = Int(hash((re_params.random_seed, i)) % typemax(Int))  # Int64, non-negative, 0 ≤ seed < 2^63  # TODO: accept any integer (Int, UInt, BigInt)?
        _, cfgs, _ = nvt_monte_carlo(configs[i], h, T, re_params.equilibrium_steps, seed)  # TODO: reproducible?
        configs[i] = cfgs[end]
        energies[i] = interacting_energy(configs[i], h).val
    end
    println("Equilibrated energies: ", energies)

    # ───────────────────────── Production ─────────────────────────
    attempted_swaps = 0
    accepted_swaps = 0

    for _ in 1:nsteps
        # 1. Local moves in every replica
        for (i, T) in enumerate(re_params.temperatures)
            cid = config_ids[i]
            _, cfgs, _ = nvt_monte_carlo(configs[i], h, T, swapint, re_params.random_seed)
            configs[i] = cfgs[end]
            energies[i] = interacting_energy(configs[i], h).val

            if store_trajectories
                append!(energy_trajs[i], [interacting_energy(c, h).val for c in cfgs])
                append!(config_trajs[i], cfgs)
                # append!(energy_trajs[cid], [interacting_energy(c, h).val for c in cfgs])
                # append!(config_trajs[cid], cfgs)
                append!(assignment_trajs[cid], fill(i, length(cfgs)))
            end
        end

        # 2. Attempt exchanges between neighbouring replicas
        for i in 1:nrep - 1
            kb = 8.617_333_262e-5  # eV K-1  # TODO: use Constants.jl
            βᵢ = 1 / (kb * re_params.temperatures[i])
            βⱼ = 1 / (kb * re_params.temperatures[i + 1])
            Δ = (βᵢ - βⱼ) * (energies[i + 1] - energies[i])

            attempted_swaps += 1
            if Δ ≤ 0 || rand() < exp(-Δ)
                energies[i], energies[i + 1] = energies[i + 1], energies[i]
                configs[i], configs[i + 1] = configs[i + 1], configs[i]
                config_ids[i], config_ids[i + 1] = config_ids[i + 1], config_ids[i]
                accepted_swaps += 1
                @debug "Swap accepted between replicas $i and $(i + 1)"
            end
        end
    end

    swap_acc_rate = accepted_swaps / max(attempted_swaps, 1)

    if store_trajectories
        return energy_trajs, config_trajs, assignment_trajs, swap_acc_rate
    else
        return energies, configs, swap_acc_rate
    end
end
