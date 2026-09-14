"""
    MuVTMCParameters <: SamplingParameters

Parameters for fixed-temperature grand-canonical (μVT) Metropolis sampling of
continuous atomistic walkers through `MC_muVT_walk!`.

# Fields
- `temperatures::Vector{Float64}`: The temperature ladder, in Kelvin.
- `activity_volumes::Vector{Float64}`: The dimensionless activity-volume zV per
  temperature (`zV = e^(βμ) V / Λ(T)^3`, folded by the caller exactly as the
  kernel documents; the driver never sees μ or Λ). Must match `temperatures` in
  length: at fixed μ the activity is temperature-dependent, so each rung carries
  its own value.
- `equilibrium_steps::Int64`: Kernel steps of equilibration per temperature, run
  in ten equal adaptation blocks (documented contract: step-size adjustment acts
  between blocks on the displacement-only acceptance rate of the block).
- `sampling_steps::Int64`: Kernel steps of production per temperature, with the
  step size frozen (adapting during production breaks detailed balance).
- `sampling_interval::Int64`: Kernel steps between recorded samples.
- `step_size::Float64`: Displacement step size (Angstrom; mutable runtime state).
- `step_size_lo::Float64`: Lower bound for the step-size adjustment.
- `step_size_up::Float64`: Upper bound for the step-size adjustment.
- `accept_range::Tuple{Float64,Float64}`: Acceptance window for the adjustment.
- `random_seed::Int64`: Base seed. Temperature rung `i` seeds equilibration with
  `random_seed + 2(i-1)` and production with `random_seed + 2(i-1) + 1`, so no
  two phases share a random stream.
"""
mutable struct MuVTMCParameters <: SamplingParameters
    temperatures::Vector{Float64}
    activity_volumes::Vector{Float64}
    equilibrium_steps::Int64
    sampling_steps::Int64
    sampling_interval::Int64
    step_size::Float64
    step_size_lo::Float64
    step_size_up::Float64
    accept_range::Tuple{Float64,Float64}
    random_seed::Int64
    function MuVTMCParameters(temperatures, activity_volumes;
                              equilibrium_steps::Int64=10_000,
                              sampling_steps::Int64=50_000,
                              sampling_interval::Int64=10,
                              step_size::Float64=0.5,
                              step_size_lo::Float64=0.01,
                              step_size_up::Float64=2.0,
                              accept_range::Tuple{Float64,Float64}=(0.25, 0.75),
                              random_seed::Int64=1234)
        if length(temperatures) != length(activity_volumes)
            throw(ArgumentError("temperatures and activity_volumes must have the same length"))
        end
        if any(t -> t <= 0.0, temperatures) || any(z -> z <= 0.0, activity_volumes)
            throw(ArgumentError("temperatures and activity_volumes must be positive"))
        end
        if sampling_interval <= 0 || equilibrium_steps < 0 || sampling_steps <= 0
            throw(ArgumentError("equilibrium_steps must be non-negative and sampling_steps, sampling_interval positive"))
        end
        new(collect(Float64, temperatures), collect(Float64, activity_volumes),
            equilibrium_steps, sampling_steps, sampling_interval,
            step_size, step_size_lo, step_size_up, accept_range, random_seed)
    end
end

"""
    monte_carlo_sampling(mc_routine::MCAtomGrandCanonicalMoves,
                         at::AtomWalker{1},
                         pot::SingleComponentPotential{Pairwise},
                         mc_params::MuVTMCParameters;
                         kb::Float64=8.617333262e-5)

Fixed-temperature grand-canonical (μVT) Metropolis sampling over a temperature
ladder, wrapping `MC_muVT_walk!`. The routine supplies the channel split
(`p_move`, `p_insert`); its nested-sampling-only fields (`step_rate_source`,
`mc_steps_per_particle`) are not consumed here. Per temperature rung:
equilibration in ten adaptation blocks (step size adjusted between blocks on the
displacement-only rate; a block that attempted no displacement skips the
adjustment), then production with the step size frozen, recording the particle
count every `sampling_interval` kernel steps and re-anchoring the walker's
incremental energy with a from-scratch `interacting_energy` recompute every ten
recorded samples. Phase seeding follows the documented per-rung law of
`MuVTMCParameters.random_seed`. The final walker of each rung starts the next
(sequential annealing over the ladder); the entry check warns once when the
potential's finite interaction range exceeds half the smallest cell edge.

# Returns
A `NamedTuple` of per-temperature results:
- `mean_N::Vector{Float64}`, `var_N::Vector{Float64}`: Particle-number moments.
- `mean_U::Vector{Float64}`: Mean total interaction energy, in eV.
- `p_N::Vector{Vector{Float64}}`, `N_support::Vector{UnitRange{Int}}`: The
  recorded particle-number histogram per rung, normalized, over `0:N_max`.
- `N_series::Vector{Vector{Int}}`: The recorded particle-number series per rung
  (the cross-check artifact; its length is `sampling_steps ÷ sampling_interval`).
- `acceptance::Vector{NamedTuple}`: Production per-channel attempt/accept
  counters per rung (the kernel's `move_stats` schema).
- `walkers::Vector{AtomWalker}`: An independent copy of the final walker per rung.
"""
function monte_carlo_sampling(mc_routine::MCAtomGrandCanonicalMoves,
                              at::AtomWalker{1},
                              pot::SingleComponentPotential{Pairwise},
                              mc_params::MuVTMCParameters;
                              kb::Float64=8.617333262e-5)
    _warn_min_image_cutoff(pot, at.configuration)
    sp = _walker_species(at)
    n_T = length(mc_params.temperatures)
    mean_N = Vector{Float64}(undef, n_T)
    var_N = Vector{Float64}(undef, n_T)
    mean_U = Vector{Float64}(undef, n_T)
    p_N = Vector{Vector{Float64}}(undef, n_T)
    N_support = Vector{UnitRange{Int}}(undef, n_T)
    N_series = Vector{Vector{Int}}(undef, n_T)
    acceptance = Vector{NamedTuple}(undef, n_T)
    walkers = Vector{typeof(at)}(undef, n_T)

    walker = deepcopy(at)
    for (i, temp) in enumerate(mc_params.temperatures)
        zV = mc_params.activity_volumes[i]
        # Equilibration: ten adaptation blocks, displacement-only rate
        Random.seed!(mc_params.random_seed + 2 * (i - 1))
        block = mc_params.equilibrium_steps ÷ 10
        for _ in 1:10
            block == 0 && break
            _, _, stats = MC_muVT_walk!(block, walker, pot, temp;
                                        zV=zV, species=sp,
                                        p_move=mc_routine.p_move,
                                        p_insert=mc_routine.p_insert,
                                        step_size=mc_params.step_size, kb=kb)
            if stats.move_attempted > 0
                adjust_step_size(mc_params, stats.move_accepted / stats.move_attempted;
                                 range=mc_params.accept_range)
            end
        end
        # Production: step size frozen, recording every sampling_interval steps
        Random.seed!(mc_params.random_seed + 2 * (i - 1) + 1)
        n_rec = mc_params.sampling_steps ÷ mc_params.sampling_interval
        ns = Vector{Int}(undef, n_rec)
        e_sum = 0.0
        att = Dict{Symbol,Int}()
        for k in 1:n_rec
            _, _, stats = MC_muVT_walk!(mc_params.sampling_interval, walker, pot, temp;
                                        zV=zV, species=sp,
                                        p_move=mc_routine.p_move,
                                        p_insert=mc_routine.p_insert,
                                        step_size=mc_params.step_size, kb=kb)
            ns[k] = walker.list_num_par[1]
            e_sum += ustrip(u"eV", walker.energy)
            for (key, val) in pairs(stats)
                att[key] = get(att, key, 0) + val
            end
            if k % 10 == 0
                # Re-anchor the incremental energy from scratch
                walker.energy = interacting_energy(walker.configuration, pot,
                    walker.list_num_par, walker.frozen) + walker.energy_frozen_part
            end
        end
        mean_N[i] = mean(ns)
        var_N[i] = var(ns)
        mean_U[i] = e_sum / n_rec
        n_max = maximum(ns)
        hist = zeros(Float64, n_max + 1)
        for n in ns
            hist[n + 1] += 1.0
        end
        p_N[i] = hist ./ n_rec
        N_support[i] = 0:n_max
        N_series[i] = ns
        acceptance[i] = (move_attempted=get(att, :move_attempted, 0),
                        move_accepted=get(att, :move_accepted, 0),
                        insert_attempted=get(att, :insert_attempted, 0),
                        insert_accepted=get(att, :insert_accepted, 0),
                        delete_attempted=get(att, :delete_attempted, 0),
                        delete_accepted=get(att, :delete_accepted, 0))
        walkers[i] = deepcopy(walker)
        @info "muVT MC T = $temp K, zV = $zV: <N> = $(round(mean_N[i]; sigdigits=5)), var(N) = $(round(var_N[i]; sigdigits=5)), <U> = $(round(mean_U[i]; sigdigits=5)) eV, step size = $(round(mc_params.step_size; sigdigits=4))"
    end

    return (mean_N=mean_N, var_N=var_N, mean_U=mean_U, p_N=p_N,
            N_support=N_support, N_series=N_series, acceptance=acceptance,
            walkers=walkers)
end

# Chemical identity for insertions, resolved once at driver entry: an empty
# single-component walker carries no species record, so sampling must start
# from a non-empty walker (mid-run empty states are fine; the identity is
# already fixed by then)
function _walker_species(walker::AtomWalker{1})
    if length(walker.configuration) == 0
        throw(ArgumentError("muVT sampling requires a non-empty starting walker (an empty single-component walker carries no species record)"))
    end
    return species(walker.configuration, 1)
end

# Catch-all for the common mistake of omitting the routine, mirroring the
# Metropolis guard (whose params annotation cannot fire for this params type)
function monte_carlo_sampling(system, potential_or_hamiltonian, mc_params::MuVTMCParameters; kwargs...)
    throw(ArgumentError(
        "`monte_carlo_sampling` requires an `MCRoutine` as the first argument.\n" *
        "For muVT sampling, use `MCAtomGrandCanonicalMoves()`.\n" *
        "Example: monte_carlo_sampling(MCAtomGrandCanonicalMoves(), walker, potential, params)"
    ))
end
