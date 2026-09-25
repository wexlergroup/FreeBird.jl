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


# ======================================================================
# Fixed-site lattice-gas μVT Metropolis sampling
# ======================================================================

function _validate_lattice_muvt_routine(mc_routine::MCGrandCanonicalMoves)
    mc_routine.p_insert > 0.0 ||
        throw(ArgumentError("lattice μVT sampling requires p_insert > 0"))
    p_delete = 1.0 - mc_routine.p_move - mc_routine.p_insert
    p_delete > 0.0 ||
        throw(ArgumentError("lattice μVT sampling requires a nonzero deletion probability"))
    mc_routine.clusters_freq == 0 ||
        throw(ArgumentError("AtomicLattice cluster moves are not implemented for μVT Metropolis sampling; use clusters_freq=0"))
    mc_routine.p_bias == 0.0 ||
        throw(ArgumentError("biased insertions are not implemented for μVT Metropolis sampling; use p_bias=0"))
    !mc_routine.incremental ||
        throw(ArgumentError("incremental lattice energies are not implemented for μVT Metropolis sampling; use incremental=false"))
    (mc_routine.p_move == 0.0 || mc_routine.swaps_freq > 0) ||
        throw(ArgumentError("p_move > 0 requires swaps_freq > 0 when clusters_freq=0"))
    return p_delete
end

function _validate_lattice_muvt_lattice(lattice::AtomicLattice{C}) where C
    AbstractWalkers._validate_atomic_components(lattice)
    C == 1 || throw(ArgumentError(
        "fixed-temperature lattice μVT sampling currently supports a " *
        "one-species AtomicLattice; got $C species"))
    return nothing
end

_lattice_muvt_energy(lattice::AtomicLattice, h::ClassicalHamiltonian) =
    ustrip(u"eV", interacting_energy(lattice, h))

_lattice_muvt_energy(lattice::AtomicLattice, calc::PyMLPotential) =
    ustrip(u"eV", interacting_energy(lattice, calc))

_lattice_muvt_copy(lattice::AtomicLattice) = deepcopy(lattice)

"""
    _lattice_muvt_acceptance_probability(ΔΩ, β, move_type, n, n_sites,
                                          p_insert, p_delete)

Internal Metropolis-Hastings acceptance probability for a uniform-site lattice
gas proposal. Insertions include
`(p_delete/p_insert) * (n_sites-n)/(n+1)` and deletions its reciprocal; fixed-N
moves are symmetric. Computing in log space avoids overflow at low temperature.
"""
function _lattice_muvt_acceptance_probability(
    delta_omega::Real,
    beta::Real,
    move_type::Symbol,
    n::Int,
    n_sites::Int,
    p_insert::Real,
    p_delete::Real,
)
    log_q_ratio = if move_type === :insert
        log(p_delete) - log(p_insert) + log(n_sites - n) - log(n + 1)
    elseif move_type === :delete
        log(p_insert) - log(p_delete) + log(n) - log(n_sites - n + 1)
    elseif move_type === :move
        0.0
    else
        throw(ArgumentError("unknown lattice μVT move type :$move_type"))
    end
    return exp(min(0.0, -beta * delta_omega + log_q_ratio))
end

"""
    metropolis_hastings(mc_routine, current_lattice, current_energy,
                        proposed_lattice, h, beta, mu, move_type)

Internal acceptance helper for one-species lattice μVT sampling.
`current_energy` and the returned energy are bare interaction energies in eV;
μ enters only through `ΔΩ = Δ(E-μN)`. The insert/delete proposal ratio is
included explicitly.
"""
function metropolis_hastings(
    mc_routine::MCGrandCanonicalMoves,
    current_lattice::AtomicLattice,
    current_energy::Float64,
    proposed_lattice::AtomicLattice,
    h::Union{ClassicalHamiltonian,PyMLPotential},
    beta::Float64,
    mu::Float64,
    move_type::Symbol,
)
    _validate_lattice_muvt_lattice(current_lattice)
    _validate_lattice_muvt_lattice(proposed_lattice)
    proposed_energy = _lattice_muvt_energy(proposed_lattice, h)
    n = n_occupied(current_lattice)
    n_proposed = n_occupied(proposed_lattice)
    delta_omega = (proposed_energy - current_energy) - mu * (n_proposed - n)
    p_delete = 1.0 - mc_routine.p_move - mc_routine.p_insert
    acceptance = _lattice_muvt_acceptance_probability(
        delta_omega, beta, move_type, n, num_sites(current_lattice),
        mc_routine.p_insert, p_delete)

    if acceptance == 1.0 || rand() < acceptance
        return proposed_lattice, proposed_energy, 1, coverage(proposed_lattice)
    end
    return current_lattice, current_energy, 0, coverage(current_lattice)
end

"""
    μvt_monte_carlo(mc_routine::MCGrandCanonicalMoves,
                    lattice::AtomicLattice{1},
                    h::Union{ClassicalHamiltonian,PyMLPotential},
                    temperature, num_steps, random_seed;
                    kb=8.617333262e-5, μ=0.0, record_interval=1)

Run fixed-temperature grand-canonical Metropolis sampling on a one-species
`AtomicLattice`. Energies are always recorded as bare E; the acceptance test
uses Δ(E-μN). Uniform insertions and deletions carry the reverse-proposal ratio
required for detailed balance. `record_interval` thins stored configurations
and always retains the final step. The supplied chemical potential must use the
same energy reference as `h`, including any per-adsorbate reference encoded by
an atomistic calculator.

Returns `(energies, configurations, coverages, accepted_steps)` at the recorded
steps. Invalid full-lattice insertions and empty-lattice deletions are guard
skips and remain in the acceptance-rate denominator.
"""
function μvt_monte_carlo(
    mc_routine::MCGrandCanonicalMoves,
    lattice::AtomicLattice,
    h::Union{ClassicalHamiltonian,PyMLPotential},
    temperature::Real,
    num_steps::Integer,
    random_seed::Integer;
    kb::Float64=8.617_333_262e-5,
    μ::Real=0.0,
    record_interval::Integer=1,
)
    temperature > 0 || throw(ArgumentError("temperature must be positive"))
    num_steps >= 0 || throw(ArgumentError("num_steps must be non-negative"))
    record_interval > 0 || throw(ArgumentError("record_interval must be positive"))
    kb > 0.0 || throw(ArgumentError("kb must be positive"))
    isfinite(μ) || throw(ArgumentError("μ must be finite"))
    _validate_lattice_muvt_lattice(lattice)
    _validate_lattice_muvt_routine(mc_routine)

    Random.seed!(random_seed)
    beta = 1.0 / (kb * Float64(temperature))
    mu = Float64(μ)
    current_lattice = _lattice_muvt_copy(lattice)
    current_energy = _lattice_muvt_energy(current_lattice, h)
    energies = Float64[]
    configurations = Vector{typeof(lattice)}()
    coverages = Float64[]
    accepted_steps = 0
    p_delete = 1.0 - mc_routine.p_move - mc_routine.p_insert

    for step in 1:num_steps
        proposed_lattice = _lattice_muvt_copy(current_lattice)
        n = n_occupied(proposed_lattice)
        r = rand()
        move_type = :none
        valid_move = true

        if r < mc_routine.p_move
            if mc_routine.swap_mode === :occupied_empty
                if n == 0 || n == num_sites(proposed_lattice)
                    valid_move = false
                else
                    from = rand(occupied_indices(proposed_lattice))
                    to = rand(empty_indices(proposed_lattice))
                    set_occupied!(proposed_lattice, from, false)
                    set_occupied!(proposed_lattice, to, true)
                end
            else
                lattice_random_walk!(proposed_lattice)
            end
            move_type = :move
        elseif r < mc_routine.p_move + mc_routine.p_insert
            valid_move, _, _ = lattice_insert_particle!(proposed_lattice)
            move_type = :insert
        else
            p_delete > 0.0 || error("unreachable zero-probability deletion branch")
            valid_move, _, _ = lattice_delete_particle!(proposed_lattice)
            move_type = :delete
        end

        if valid_move
            current_lattice, current_energy, accepted, _ = metropolis_hastings(
                mc_routine, current_lattice, current_energy, proposed_lattice,
                h, beta, mu, move_type)
            accepted_steps += accepted
        end

        if step % record_interval == 0 || step == num_steps
            push!(energies, current_energy)
            push!(configurations, _lattice_muvt_copy(current_lattice))
            push!(coverages, coverage(current_lattice))
        end
    end

    return energies, configurations, coverages, accepted_steps
end

"""
    monte_carlo_sampling(mc_routine::MCGrandCanonicalMoves,
                         lattice::AtomicLattice{1},
                         h::Union{ClassicalHamiltonian,PyMLPotential},
                         mc_params::MetropolisMCParameters;
                         kb=8.617333262e-5, sampling_interval=1)

Run one-species lattice μVT Metropolis sampling over every chemical potential
and temperature in `mc_params`. Each μ starts from the supplied lattice; within
a μ the final configuration at one temperature seeds the next. Equilibration
and production use distinct deterministic seeds for every `(μ,T)` rung.

Returns a result `DataFrame` with bare mean energy, `c_omega`, coverage and
acceptance rate, plus a dictionary containing the final lattice for each μ.
"""
function _lattice_muvt_sampling(
    mc_routine::MCGrandCanonicalMoves,
    lattice::AtomicLattice,
    h::Union{ClassicalHamiltonian,PyMLPotential},
    mc_params::MetropolisMCParameters;
    kb::Float64=8.617_333_262e-5,
    sampling_interval::Integer=1,
)
    mus = mc_params.chemical_potentials
    isnothing(mus) && throw(ArgumentError(
        "MCGrandCanonicalMoves requires chemical_potentials in MetropolisMCParameters"))
    mc_params.equilibrium_steps >= 0 ||
        throw(ArgumentError("equilibrium_steps must be non-negative"))
    mc_params.sampling_steps > 0 ||
        throw(ArgumentError("sampling_steps must be positive"))
    sampling_interval > 0 ||
        throw(ArgumentError("sampling_interval must be positive"))
    _validate_lattice_muvt_lattice(lattice)
    _validate_lattice_muvt_routine(mc_routine)

    rows = NamedTuple[]
    final_configs = Dict{Float64,typeof(lattice)}()
    n_temperatures = length(mc_params.temperatures)

    for (mu_index, mu) in enumerate(mus)
        current_lattice = _lattice_muvt_copy(lattice)
        for (temp_index, temp) in enumerate(mc_params.temperatures)
            rung = (mu_index - 1) * n_temperatures + temp_index
            eq_seed = mc_params.random_seed + 2 * (rung - 1)
            if mc_params.equilibrium_steps > 0
                _, eq_configs, _, _ = μvt_monte_carlo(
                    mc_routine, current_lattice, h, temp,
                    mc_params.equilibrium_steps, eq_seed;
                    kb=kb, μ=mu,
                    record_interval=max(mc_params.equilibrium_steps, 1))
                current_lattice = eq_configs[end]
            end

            sample_energies, sample_configs, sample_coverages, accepted =
                μvt_monte_carlo(
                    mc_routine, current_lattice, h, temp,
                    mc_params.sampling_steps, eq_seed + 1;
                    kb=kb, μ=mu, record_interval=sampling_interval)
            particle_counts = n_occupied.(sample_configs)
            omega = sample_energies .- mu .* particle_counts
            c_omega = var(omega; corrected=false) / (kb * temp^2)
            acceptance_rate = accepted / mc_params.sampling_steps
            current_lattice = sample_configs[end]

            push!(rows, (
                temperature=Float64(temp),
                energy=mean(sample_energies),
                c_omega=c_omega,
                cov=mean(sample_coverages),
                acceptance_rate=acceptance_rate,
                chemical_potential=mu,
            ))
        end
        final_configs[mu] = _lattice_muvt_copy(current_lattice)
    end

    return DataFrame(rows), final_configs
end

function monte_carlo_sampling(
    mc_routine::MCGrandCanonicalMoves,
    lattice::AtomicLattice,
    h::ClassicalHamiltonian,
    mc_params::MetropolisMCParameters;
    kb::Float64=8.617_333_262e-5,
    sampling_interval::Integer=1,
)
    return _lattice_muvt_sampling(
        mc_routine, lattice, h, mc_params;
        kb=kb, sampling_interval=sampling_interval)
end

function monte_carlo_sampling(
    mc_routine::MCGrandCanonicalMoves,
    lattice::AtomicLattice,
    calc::PyMLPotential,
    mc_params::MetropolisMCParameters;
    kb::Float64=8.617_333_262e-5,
    sampling_interval::Integer=1,
)
    return _lattice_muvt_sampling(
        mc_routine, lattice, calc, mc_params;
        kb=kb, sampling_interval=sampling_interval)
end
