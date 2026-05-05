# ======================================================================
# Grand-canonical nested sampling — U-sorted, ideal-gas-reference variant
#
# This file implements a second variant of grand-canonical nested sampling
# for lattice systems. Distinct from the Ω-sorted construction in
# `nested_sampling.jl`, here the chemical potential μ enters the prior via
# a reference fugacity z₀ rather than the sort quantity. A single run
# covers a neighborhood of μ that can be reweighted in post-processing.
#
# Algorithmic core:
#   - Prior:    π₀(σ) ∝ z₀^N(σ); per-site Bernoulli with p₀ = z₀/(1+z₀).
#   - Sort:     U(σ) — the potential energy. μ does not enter the loop.
#   - Moves:    single-site flip (Metropolis ratio min(1, z₀^ΔN))
#               + fixed-N moves (local swap + geometric cluster).
#
# Both constructions are intended to coexist as parallel sampling schemes;
# this file does not modify the Ω-sorted code.
# ======================================================================

"""
    struct MCIdealGasReferenceMoves <: MCRoutine

A type for generating a new walker using U-sorted, ideal-gas-reference
grand-canonical MCMC moves: single-site flips with Metropolis ratio
``\\min(1, z_0^{\\Delta N})`` against the prior ``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}``,
mixed with fixed-N moves (local swaps and/or geometric cluster moves).

The chemical potential ``\\mu`` is *not* a sampling parameter — it enters
only during post-processing reweighting via importance weights
``(z/z_0)^{N_j}`` where ``z = e^{\\beta \\mu}``.

# Fields
- `p_flip::Float64`: Probability of a single-site flip per MCMC step (default 0.5).
  The remaining `1 - p_flip` is the probability of a fixed-N move.
- `swaps_freq::Int`: Relative weight of local swaps within the fixed-N branch (default 1).
- `clusters_freq::Int`: Relative weight of cluster moves within the fixed-N branch (default 0 = disabled).
- `initial_cluster_p::Float64`: Starting growth probability for geometric cluster moves (default 0.3).
- `target_cluster_accept::Float64`: Target acceptance rate for adaptive cluster p tuning (default 0.3).
- `cluster_adjust_interval::Int`: NS iterations between cluster p adjustments (default 50).
- `cluster_p_floor::Float64`: Lower bound for adaptive cluster p (default 0.01).
- `cluster_p_ceiling::Float64`: Upper bound for adaptive cluster p (default 1.0).

When `clusters_freq == 0` (default), the fixed-N branch performs only local
swaps. When `clusters_freq > 0`, fixed-N moves are mixed by the
`swaps_freq:clusters_freq` ratio.
"""
struct MCIdealGasReferenceMoves <: MCRoutine
    p_flip::Float64
    swaps_freq::Int
    clusters_freq::Int
    initial_cluster_p::Float64
    target_cluster_accept::Float64
    cluster_adjust_interval::Int
    cluster_p_floor::Float64
    cluster_p_ceiling::Float64
    function MCIdealGasReferenceMoves(;
            p_flip::Float64=0.5,
            swaps_freq::Int=1,
            clusters_freq::Int=0,
            initial_cluster_p::Float64=0.3,
            target_cluster_accept::Float64=0.3,
            cluster_adjust_interval::Int=50,
            cluster_p_floor::Float64=0.01,
            cluster_p_ceiling::Float64=1.0)
        if !(0.0 <= p_flip <= 1.0)
            throw(ArgumentError("p_flip must be in [0,1], got $p_flip"))
        end
        if swaps_freq < 0 || clusters_freq < 0
            throw(ArgumentError("swaps_freq and clusters_freq must be non-negative"))
        end
        if p_flip < 1.0 && swaps_freq + clusters_freq <= 0
            throw(ArgumentError("If p_flip < 1, swaps_freq + clusters_freq must be > 0"))
        end
        if !(0.0 <= initial_cluster_p <= 1.0)
            throw(ArgumentError("initial_cluster_p must be in [0,1], got $initial_cluster_p"))
        end
        if !(0.0 <= cluster_p_floor <= cluster_p_ceiling <= 1.0)
            throw(ArgumentError("require 0 <= cluster_p_floor <= cluster_p_ceiling <= 1"))
        end
        new(p_flip, swaps_freq, clusters_freq,
            initial_cluster_p, target_cluster_accept, cluster_adjust_interval,
            cluster_p_floor, cluster_p_ceiling)
    end
end

"""
    mutable struct LatticeGCNSIdealGasReferenceParameters <: SamplingParameters

Parameters for U-sorted, ideal-gas-reference grand-canonical nested sampling
on lattice systems.

Walkers are initialized from the prior ``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}``
(per-site Bernoulli with ``p_0 = z_0/(1+z_0)``) and the sort quantity is
the potential energy ``U(\\sigma)``. The chemical potential ``\\mu`` does
**not** appear in the sampling loop — it is supplied to post-processing
functions to reweight the run output.

# Fields
- `mc_steps::Int64`: MCMC steps per replacement walker.
- `z0::Float64`: Reference fugacity (dimensionless), `z0 > 0`. Defines the
  prior ``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}``. The case `z0 = 1`
  recovers the uniform prior over microstates (matching the prior implicit
  in the existing Ω-sorted GCNS).
- `energy_perturbation::Float64`: Perturbation magnitude for breaking energy degeneracies.
- `random_seed::Int64`: Seed for the random number generator.
- `fail_count::Int64`: Consecutive failed replacements (mutable runtime state).
- `allowed_fail_count::Int64`: Maximum consecutive failures before a warning is emitted.
- `n_max::Int64`: Upper bound on particle count per walker. Default `typemax(Int64)`.
- `T0_ref::Union{Float64,Nothing}`: Optional reference temperature (K). When set,
  used as the default `T0` argument by [`effective_sample_size_ideal_gas_reference`](@ref);
  has no effect on the sampling loop. Default `nothing`.
- `cluster_p::Float64`: Current cluster growth probability (mutable runtime state).
- `cluster_accepted::Float64`: Accepted cluster moves in the current adjustment window.
- `cluster_total::Float64`: Total cluster moves attempted in the current adjustment window.
- `cluster_p_history::Vector{Float64}`: Trajectory of `cluster_p` after each adjustment.
- `cluster_accept_history::Vector{Float64}`: Acceptance rate at each adjustment.
- `cluster_adjust_iterations::Vector{Int}`: NS iteration index at each adjustment.
"""
mutable struct LatticeGCNSIdealGasReferenceParameters <: SamplingParameters
    mc_steps::Int64
    z0::Float64
    energy_perturbation::Float64
    random_seed::Int64
    fail_count::Int64
    allowed_fail_count::Int64
    n_max::Int64
    T0_ref::Union{Float64,Nothing}
    cluster_p::Float64
    cluster_accepted::Float64
    cluster_total::Float64
    cluster_p_history::Vector{Float64}
    cluster_accept_history::Vector{Float64}
    cluster_adjust_iterations::Vector{Int}
end

"""
    LatticeGCNSIdealGasReferenceParameters(; mc_steps=100, z0=1.0, ...)

Convenience keyword constructor for [`LatticeGCNSIdealGasReferenceParameters`](@ref).
See the type docstring for the full field list and defaults.
"""
function LatticeGCNSIdealGasReferenceParameters(;
        mc_steps::Int64=100,
        z0::Float64=1.0,
        energy_perturbation::Float64=1e-12,
        random_seed::Int64=1234,
        fail_count::Int64=0,
        allowed_fail_count::Int64=10,
        n_max::Int64=typemax(Int64),
        T0_ref::Union{Float64,Nothing}=nothing,
        cluster_p::Float64=0.3,
        cluster_accepted::Float64=0.0,
        cluster_total::Float64=0.0,
        cluster_p_history::Vector{Float64}=Float64[],
        cluster_accept_history::Vector{Float64}=Float64[],
        cluster_adjust_iterations::Vector{Int}=Int[])
    if z0 <= 0.0
        throw(ArgumentError("z0 must be positive, got $z0"))
    end
    if T0_ref !== nothing && T0_ref <= 0.0
        throw(ArgumentError("T0_ref must be positive, got $T0_ref"))
    end
    LatticeGCNSIdealGasReferenceParameters(
        mc_steps, z0, energy_perturbation, random_seed,
        fail_count, allowed_fail_count, n_max, T0_ref,
        cluster_p, cluster_accepted, cluster_total,
        cluster_p_history, cluster_accept_history, cluster_adjust_iterations)
end

# ======================================================================
# Walker initialization
# ======================================================================

"""
    _init_ideal_gas_reference_walkers!(liveset::LatticeGasWalkers,
                                        params::LatticeGCNSIdealGasReferenceParameters)

Initialize each walker's configuration by drawing from the ideal-gas
reference prior ``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}``. Each site is
occupied independently with probability ``p_0 = z_0 / (1 + z_0)``.

If a walker exceeds `params.n_max` occupied sites after the Bernoulli draw,
particles are randomly deleted until the cap is satisfied.

Each walker's energy is reassigned via [`assign_energy!`](@ref) using the
liveset's Hamiltonian and `params.energy_perturbation`.
"""
function _init_ideal_gas_reference_walkers!(liveset::LatticeGasWalkers,
                                             params::LatticeGCNSIdealGasReferenceParameters)
    h = liveset.hamiltonian
    p0 = params.z0 / (1.0 + params.z0)
    n_max = params.n_max
    for walker in liveset.walkers
        random_microstate!(walker.configuration; p=p0)
        n_occ = sum(walker.configuration.components[1])
        if n_occ > n_max
            occupied = findall(walker.configuration.components[1])
            shuffle!(occupied)
            for i in 1:(n_occ - n_max)
                walker.configuration.components[1][occupied[i]] = false
            end
        end
        assign_energy!(walker, h; perturb_energy=params.energy_perturbation)
    end
    return liveset
end

# ======================================================================
# Sampling step and top-level loop
# ======================================================================

"""
    nested_sampling_step!(liveset::LatticeGasWalkers,
                          params::LatticeGCNSIdealGasReferenceParameters,
                          mc_routine::MCIdealGasReferenceMoves;
                          ns_iteration::Int=0)

Perform one step of U-sorted, ideal-gas-reference grand-canonical nested
sampling.

Sorts walkers by potential energy ``U`` (descending), removes the worst
(highest ``U``), clones a parent with ``U < U_{\\max}``, and decorrelates
the clone via [`MC_ideal_gas_reference_walk!`](@ref). Cluster-move
acceptance counts are accumulated and `params.cluster_p` is adapted via
[`adjust_cluster_p`](@ref) every `mc_routine.cluster_adjust_interval`
steps.

# Returns
- `iter::Union{Missing,Int}`: NS iteration number, or `missing` if the step failed.
- `u_worst`: The ``U`` value of the removed walker (with units).
- `n_worst::Union{Missing,Int}`: The particle count of the removed walker.
- `liveset`: The updated liveset.
- `params`: The updated parameters.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               params::LatticeGCNSIdealGasReferenceParameters,
                               mc_routine::MCIdealGasReferenceMoves;
                               ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    n_walkers = length(ats)

    iter::Union{Missing,Int} = missing
    worst = ats[1]
    u_worst = worst.energy
    n_worst::Union{Missing,Int} = sum(worst.configuration.components[1])

    u_max_val = u_worst.val

    # Prefer parents strictly below u_worst to avoid stalling on tied energies
    eligible = [k for k in 2:n_walkers if ats[k].energy < u_worst]
    parent_idx = isempty(eligible) ? rand(2:n_walkers) : rand(eligible)
    to_walk = deepcopy(ats[parent_idx])

    accept, _, to_walk, cl_accepted, cl_total = MC_ideal_gas_reference_walk!(
        params.mc_steps, to_walk, h, u_max_val, params.z0;
        p_flip=mc_routine.p_flip,
        swaps_freq=mc_routine.swaps_freq,
        clusters_freq=mc_routine.clusters_freq,
        cluster_p=params.cluster_p,
        n_max=params.n_max,
        energy_perturb=params.energy_perturbation)

    if accept
        push!(ats, to_walk)
        popfirst!(ats)
        update_iter!(liveset)
        params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        u_worst = missing
        n_worst = missing
        params.fail_count += 1
    end

    if cl_total > 0
        params.cluster_accepted += cl_accepted
        params.cluster_total += cl_total
        if mc_routine.cluster_adjust_interval > 0 &&
           params.cluster_total >= mc_routine.cluster_adjust_interval
            window_rate = params.cluster_accepted / max(params.cluster_total, 1.0)
            adjust_cluster_p(params, window_rate, ns_iteration;
                             target=mc_routine.target_cluster_accept,
                             floor=mc_routine.cluster_p_floor,
                             ceiling=mc_routine.cluster_p_ceiling)
            params.cluster_accepted = 0.0
            params.cluster_total = 0.0
        end
    end

    return iter, u_worst, n_worst, liveset, params
end

"""
    ideal_gas_reference_grand_canonical_nested_sampling(
        liveset::LatticeGasWalkers,
        params::LatticeGCNSIdealGasReferenceParameters,
        n_steps::Int64,
        mc_routine::MCIdealGasReferenceMoves,
        save_strategy::DataSavingStrategy)

Run the U-sorted, ideal-gas-reference grand-canonical nested sampling loop.

Walkers are first re-initialized from the prior
``\\pi_0(\\sigma) \\propto z_0^{N(\\sigma)}`` (per-site Bernoulli with
``p_0 = z_0/(1+z_0)``); any incoming walker state is overwritten. The loop
then iterates: sort by ``U``, remove the worst walker, record
``(\\mathrm{iter}, U, N)``, replace it with a decorrelated clone of a
surviving walker.

Post-processing of the returned DataFrame (e.g.
[`ideal_gas_reference_thermodynamic_stats`](@ref)) reweights the run to
arbitrary ``(\\mu, T)`` via importance weights ``(z/z_0)^N``.

# Arguments
- `liveset::LatticeGasWalkers`: The initial liveset (walkers will be re-initialized from the prior).
- `params::LatticeGCNSIdealGasReferenceParameters`: Sampling parameters including the reference fugacity `z0`.
- `n_steps::Int64`: Number of NS iterations.
- `mc_routine::MCIdealGasReferenceMoves`: The MCMC move recipe.
- `save_strategy::DataSavingStrategy`: Periodic-output strategy.

# Returns
- `df::DataFrame`: Columns `[:iter, :emax, :num_particles]`. `emax` is the
  potential energy of the culled walker at each successful iteration (in
  the Hamiltonian's energy units, stored as a unitless `Float64`).
- `liveset::LatticeGasWalkers`: The final surviving walkers.
- `params::LatticeGCNSIdealGasReferenceParameters`: Updated parameters
  (including final `cluster_p` and adaptation history).
"""
function ideal_gas_reference_grand_canonical_nested_sampling(
        liveset::LatticeGasWalkers,
        params::LatticeGCNSIdealGasReferenceParameters,
        n_steps::Int64,
        mc_routine::MCIdealGasReferenceMoves,
        save_strategy::DataSavingStrategy)
    _init_ideal_gas_reference_walkers!(liveset, params)

    if mc_routine.clusters_freq > 0
        params.cluster_p = mc_routine.initial_cluster_p
        params.cluster_accepted = 0.0
        params.cluster_total = 0.0
        empty!(params.cluster_p_history)
        empty!(params.cluster_accept_history)
        empty!(params.cluster_adjust_iterations)
    end

    df = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[])

    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)

        iter, u_worst, n_worst, liveset, params = nested_sampling_step!(
            liveset, params, mc_routine; ns_iteration=i)

        @debug "Ideal-gas-ref GCNS step $i, iter: $iter, U: $u_worst, N: $n_worst"

        if params.fail_count >= params.allowed_fail_count
            @warn "Ideal-gas-ref GCNS: Failed $(params.allowed_fail_count) times in a row."
            params.fail_count = 0
        end

        if !(iter isa typeof(missing))
            push!(df, (iter, u_worst.val, n_worst))
        end

        if print_info && !(iter isa typeof(missing))
            @info "Ideal-gas-ref GCNS iter: $(iter), U: $(u_worst), N: $(n_worst)"
        elseif print_info && iter isa typeof(missing)
            @info "Ideal-gas-ref GCNS MC move failed, step: $(i)"
        end

        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end

    return df, liveset, params
end
