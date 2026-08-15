"""
    mutable struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in the nested sampling scheme.

# Fields
- `mc_steps::Int64`: The number of total Monte Carlo moves to perform. For a parallel MC routine, this number will be distributed among workers.
If `mc_steps` is not divisible by the number of workers, the actual number of MC moves per worker will be `ceil(mc_steps / nworkers())`.
- `initial_step_size::Float64`: The initial step size, which is the fallback step size if MC routine fails to accept a move.
- `step_size::Float64`: The on-the-fly step size used in the sampling process.
- `step_size_lo::Float64`: The lower bound of the step size.
- `step_size_up::Float64`: The upper bound of the step size.
- `accept_range::Tuple{Float64, Float64}`: The range of acceptance rates for adjusting the step size.
e.g. (0.25, 0.75) means that the step size will decrease if the acceptance rate is below 0.25 and increase if it is above 0.75.
- `fail_count::Int64`: The number of failed MC moves in a row.
- `allowed_fail_count::Int64`: The maximum number of failed MC moves allowed before resetting the step size.
- `energy_perturbation::Float64`: The perturbation value used to adjust the energy of the walkers.
- `random_seed::Int64`: The seed for the random number generator.
- `cluster_p::Float64`: Current cluster growth probability for geometric cluster moves (mutable runtime state).
- `cluster_accepted::Float64`: Accepted cluster moves in the current adjustment window.
- `cluster_total::Float64`: Total cluster moves attempted in the current adjustment window.
- `cluster_p_history::Vector{Float64}`: Trajectory of cluster_p values after each adaptive adjustment.
- `cluster_accept_history::Vector{Float64}`: Acceptance rate at each adaptive adjustment.
- `cluster_adjust_iterations::Vector{Int}`: NS iteration index at each adaptive adjustment.
- `plateau_refill_target::Int64`: Internal bookkeeping for plateau-aware culling in the
  serial atomistic steps: the live-set size to restore once a block of exact energy ties
  has been fully evicted without replacement. `0` (the default) means no plateau block is
  in progress; the field is set and cleared by `nested_sampling_step!` and should not be
  set by callers.
"""
mutable struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    initial_step_size::Float64
    step_size::Float64
    step_size_lo::Float64
    step_size_up::Float64
    accept_range::Tuple{Float64, Float64}
    fail_count::Int64
    allowed_fail_count::Int64
    energy_perturbation::Float64
    random_seed::Int64
    cluster_p::Float64
    cluster_accepted::Float64
    cluster_total::Float64
    cluster_p_history::Vector{Float64}
    cluster_accept_history::Vector{Float64}
    cluster_adjust_iterations::Vector{Int}
    plateau_refill_target::Int64
end

function NestedSamplingParameters(;
            mc_steps::Int64=200,
            initial_step_size::Float64=0.01,
            step_size::Float64=0.1,
            step_size_lo::Float64=1e-6,
            step_size_up::Float64=1.0,
            accept_range::Tuple{Float64, Float64}=(0.25, 0.75),
            fail_count::Int64=0,
            allowed_fail_count::Int64=100,
            energy_perturbation::Float64=1e-12,
            random_seed::Int64=1234,
            cluster_p::Float64=0.3,
            cluster_accepted::Float64=0.0,
            cluster_total::Float64=0.0,
            cluster_p_history::Vector{Float64}=Float64[],
            cluster_accept_history::Vector{Float64}=Float64[],
            cluster_adjust_iterations::Vector{Int}=Int[],
            plateau_refill_target::Int64=0,
            )
    NestedSamplingParameters(mc_steps, initial_step_size, step_size, step_size_lo, step_size_up, accept_range, fail_count, allowed_fail_count, energy_perturbation, random_seed, cluster_p, cluster_accepted, cluster_total, cluster_p_history, cluster_accept_history, cluster_adjust_iterations, plateau_refill_target)
end

"""
    LatticeNestedSamplingParameters(;
            mc_steps::Int64=100,
            energy_perturbation::Float64=1e-12,
            fail_count::Int64=0,
            allowed_fail_count::Int64=10,
            random_seed::Int64=1234,
            cluster_p::Float64=0.3,
            )
A convenience constructor for `NestedSamplingParameters` with default values suitable for lattice systems.
"""
function LatticeNestedSamplingParameters(;
            mc_steps::Int64=100,
            energy_perturbation::Float64=1e-12,
            fail_count::Int64=0,
            allowed_fail_count::Int64=10,
            random_seed::Int64=1234,
            cluster_p::Float64=0.3,
            )
    NestedSamplingParameters(mc_steps=mc_steps, fail_count=fail_count, allowed_fail_count=allowed_fail_count, energy_perturbation=energy_perturbation, random_seed=random_seed, cluster_p=cluster_p)
end


"""
    abstract type MCRoutine

An abstract type representing a Monte Carlo routine.

Currently, the following concrete types are supported:
- `MCRandomWalkMaxE`: A type for generating a new walker by performing a random walk for decorrelation on the
highest-energy walker.
- `MCRandomWalkClone`: A type for generating a new walker by cloning an existing walker and performing a random walk
for decorrelation.
- `MCNewSample`: A type for generating a new walker from a random configuration. Currently, it is intended to use 
this routine for lattice gas systems.
- `MCMixedMoves`: A type for generating a new walker by performing random walks and swapping atoms. Currently, it is
intended to use this routine for multi-component systems. The actual number of random walks and swaps to perform is
determined by the weights of the fields `walks_freq` and `swaps_freq`. See [`MCMixedMoves`](@ref).
- `MCRejectionSampling`: A type for generating a new walker by performing rejection sampling. Currently, it is intended
to use this routine for lattice gas systems.
- `MCDistributed`: A type for generating new walkers by performing random walks for decorrelation in parallel using Distributed.jl.
This routine supports multiple culling walkers and multiple decorrelation walkers. See [`MCDistributed`](@ref).
"""
abstract type MCRoutine end

"""
    abstract type MCRoutineParallel <: MCRoutine
(Internal) An abstract type representing a parallel Monte Carlo routine.
"""
abstract type MCRoutineParallel <: MCRoutine end

"""
    struct MCRandomWalkMaxE <: MCRoutine
A type for generating a new walker by performing a random walk for decorrelation on the highest-energy walker.
"""
struct MCRandomWalkMaxE <: MCRoutine 
    dims::Vector{Int64}
    function MCRandomWalkMaxE(dims::Vector{Int64}=[1, 2, 3])
        new(dims)
    end
end

"""
    struct MCRandomWalkClone <: MCRoutine
A type for generating a new walker by cloning an existing walker and performing a random walk for decorrelation.
"""
struct MCRandomWalkClone <: MCRoutine 
    dims::Vector{Int64}
    function MCRandomWalkClone(;dims::Vector{Int64}=[1, 2, 3])
        new(dims)
    end
end

"""
    struct MCRandomWalkCloneParallel <: MCRoutineParallel
A type for generating a new walker by cloning an existing walker and performing a random walk for decorrelation in parallel.
"""
struct MCRandomWalkCloneParallel <: MCRoutineParallel
    dims::Vector{Int64}
    function MCRandomWalkCloneParallel(;dims::Vector{Int64}=[1, 2, 3])
        new(dims)
    end
end

"""
    struct MCDistributed <: MCRoutineParallel
A type for generating new walkers by performing random walks for decorrelation in parallel using Distributed.jl.
# Fields
- `n_cull::Int64`: The number of lowest-energy walkers to cull (replace) in each iteration. The default is 1.
- `n_decorr::Int64`: The number of walkers to use for decorrelation (random walks). The default is `nworkers() - 1`.
- `dims::Vector{Int64}`: The dimensions along which to perform the random walks.
"""
struct MCDistributed <: MCRoutineParallel
    n_cull::Int64
    n_decorr::Int64
    dims::Vector{Int64}
    function MCDistributed(;n_cull::Int64=1, n_decorr::Int64=nworkers()-1, dims::Vector{Int64}=[1, 2, 3])
        if n_cull + n_decorr != nworkers()
            error("n_cull + n_decorr must be equal to the number of workers: $(nworkers())")
        end
        @info "Distributed nested sampling initiated: n_cull: $n_cull, n_decorr: $n_decorr, total workers: $(n_cull + n_decorr)"
        new(n_cull, n_decorr, dims)
    end
end

"""
    MCRandomWalkMaxEParallel <: MCRoutineParallel
A type for generating a new walker by performing a random walk for decorrelation on the highest-energy walker(s) in parallel.
"""
struct MCRandomWalkMaxEParallel <: MCRoutineParallel
    dims::Vector{Int64}
    function MCRandomWalkMaxEParallel(;dims::Vector{Int64}=[1, 2, 3])
        new(dims)
    end
end

"""
    struct MCNewSample <: MCRoutine
A type for generating a new walker from a random configuration. Currently, it is intended to use this routine for lattice gas systems.
"""
struct MCNewSample <: MCRoutine end

"""
    struct MCMixedMoves <: MCRoutine
A type for generating a new walker by performing a mix of random walks, atom swaps, and/or geometric cluster moves.
For atomistic systems, the `walks_freq` and `swaps_freq` fields control the ratio of random walks to atom swaps.
For lattice systems, `walks_freq` and `clusters_freq` control the ratio of local swap moves to geometric cluster moves.

# Fields
- `walks_freq::Int`: The frequency of random walks (atomistic) or local swaps (lattice) to perform.
- `swaps_freq::Int`: The frequency of atom swaps to perform (atomistic only).
- `clusters_freq::Int`: The frequency of geometric cluster moves to perform (lattice only, default 0).
- `initial_cluster_p::Float64`: Starting growth probability for cluster moves (default 0.3).
- `target_cluster_accept::Float64`: Target acceptance rate for adaptive cluster p tuning (default 0.3).
- `cluster_adjust_interval::Int`: Number of NS iterations between cluster p adjustments (default 50).
- `cluster_p_floor::Float64`: Lower bound for adaptive cluster p (default 0.01).
- `cluster_p_ceiling::Float64`: Upper bound for adaptive cluster p (default 1.0).
"""
mutable struct MCMixedMoves <: MCRoutine
    walks_freq::Int
    swaps_freq::Int
    clusters_freq::Int
    initial_cluster_p::Float64
    target_cluster_accept::Float64
    cluster_adjust_interval::Int
    cluster_p_floor::Float64
    cluster_p_ceiling::Float64
end

# Backward-compatible constructor: MCMixedMoves(5, 1)
function MCMixedMoves(walks_freq::Int, swaps_freq::Int)
    MCMixedMoves(walks_freq, swaps_freq, 0, 0.3, 0.3, 50, 0.01, 1.0)
end

# Keyword constructor for lattice use
function MCMixedMoves(;
    walks_freq::Int=1,
    swaps_freq::Int=0,
    clusters_freq::Int=1,
    initial_cluster_p::Float64=0.3,
    target_cluster_accept::Float64=0.3,
    cluster_adjust_interval::Int=50,
    cluster_p_floor::Float64=0.01,
    cluster_p_ceiling::Float64=1.0,
)
    MCMixedMoves(walks_freq, swaps_freq, clusters_freq,
                 initial_cluster_p, target_cluster_accept, cluster_adjust_interval,
                 cluster_p_floor, cluster_p_ceiling)
end

"""
    struct MCMixedMovesParallel <: MCRoutineParallel
A type for generating a new walker by performing random walks and swapping atoms in parallel. Currently, it is intended to use this routine for
multi-component systems. The actual number of random walks and swaps to perform is determined by the weights of the fields `walks_freq` and `swaps_freq`.
For example, if `walks_freq=4` and `swaps_freq=1`, then the probability of performing a random walk is 4/5, and the probability of performing a swap is 1/5.

# Fields
- `walks_freq::Int`: The frequency of random walks to perform.
- `swaps_freq::Int`: The frequency of atom swaps to perform.
"""
mutable struct MCMixedMovesParallel <: MCRoutineParallel
    walks_freq::Int
    swaps_freq::Int
end

"""
    struct MCRejectionSampling <: MCRoutine
A type for generating a new walker by performing rejection sampling. Currently, it is intended to use this routine for lattice gas systems.
"""
struct MCRejectionSampling <: MCRoutine end


# ======================================================================
# Grand-canonical nested sampling types
# ======================================================================

"""
    struct MCGrandCanonicalMoves <: MCRoutine

A type for generating a new walker using grand-canonical MCMC moves that mix
fixed-N moves (local swaps and/or geometric cluster moves) with single-site
particle insertion and deletion.

# Fields
- `p_move::Float64`: Probability of a fixed-N move per MCMC step (default 0.5).
- `p_insert::Float64`: Probability of a particle insertion per step (default 0.25).
  The deletion probability is `1 - p_move - p_insert`.
- `clusters_freq::Int`: Relative weight of cluster moves within the fixed-N branch (default 0 = disabled).
- `swaps_freq::Int`: Relative weight of local swaps within the fixed-N branch (default 1).
- `initial_cluster_p::Float64`: Starting growth probability for geometric cluster moves (default 0.3).
- `target_cluster_accept::Float64`: Target acceptance rate for adaptive cluster p tuning (default 0.3).
- `cluster_adjust_interval::Int`: Number of NS iterations between cluster p adjustments (default 50).
- `cluster_p_floor::Float64`: Lower bound for adaptive cluster p (default 0.01).
- `cluster_p_ceiling::Float64`: Upper bound for adaptive cluster p (default 1.0).
- `p_bias::Float64`: Probability an insertion draws from the biased site set
  (default 0.0 = legacy uniform insertions).
- `bias_predicate::Symbol`: Biased-set predicate, `:contact` or `:cavity`
  (default `:contact`).
- `bias_shells::Int`: Neighbor shells scanned by the predicate (default 1).
- `incremental::Bool`: Opt-in incremental energy evaluation in the walk kernel
  (default `false` = the shipped full-recompute arithmetic, draw-count and
  digit identical). When `true`, non-cluster proposals under a Hamiltonian
  with `supports_site_deltas` advance a per-walk raw-energy anchor by exact
  O(z) `site_flip_delta` sums; cluster proposals and unsupported Hamiltonians
  fall back to the full recompute, which also re-anchors the accumulator.
  The delta path accumulates energy in a different floating-point order, so
  same-seed trajectories are not digit-identical to the default; flipping
  the default is deliberately out of scope (a published-run reproducibility
  contract).

When `clusters_freq == 0` (the default), the fixed-N branch uses only local swaps
(`lattice_random_walk!`), preserving backward compatibility with existing scripts.
When `clusters_freq > 0`, the fixed-N branch mixes geometric cluster moves with
local swaps according to the `clusters_freq:swaps_freq` ratio.

With `p_bias > 0` insertions mix a sub-channel restricted to
`lattice_biased_sites(x; predicate=bias_predicate, shells=bias_shells)` into
the uniform proposal; the walk's composite Metropolis-Hastings correction
keeps the sampled prior unchanged. `p_bias = 1.0` constructs but warns: the
pure biased channel freezes N whenever the biased set is empty, so the
uniform sub-channel is what repairs ergodicity.
"""
struct MCGrandCanonicalMoves <: MCRoutine
    p_move::Float64
    p_insert::Float64
    clusters_freq::Int
    swaps_freq::Int
    initial_cluster_p::Float64
    target_cluster_accept::Float64
    cluster_adjust_interval::Int
    cluster_p_floor::Float64
    cluster_p_ceiling::Float64
    p_bias::Float64
    bias_predicate::Symbol
    bias_shells::Int
    incremental::Bool
    function MCGrandCanonicalMoves(;
            p_move::Float64=0.5,
            p_insert::Float64=0.25,
            clusters_freq::Int=0,
            swaps_freq::Int=1,
            initial_cluster_p::Float64=0.3,
            target_cluster_accept::Float64=0.3,
            cluster_adjust_interval::Int=50,
            cluster_p_floor::Float64=0.01,
            cluster_p_ceiling::Float64=1.0,
            p_bias::Float64=0.0,
            bias_predicate::Symbol=:contact,
            bias_shells::Int=1,
            incremental::Bool=false)
        if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
            throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
        end
        if !(0.0 <= p_bias <= 1.0)
            throw(ArgumentError("p_bias must satisfy 0 <= p_bias <= 1"))
        end
        if bias_predicate !== :contact && bias_predicate !== :cavity
            throw(ArgumentError("unknown bias_predicate :$bias_predicate; expected :contact or :cavity"))
        end
        if bias_shells < 1
            throw(ArgumentError("bias_shells must be >= 1, got $bias_shells"))
        end
        if p_bias == 1.0
            @warn "p_bias = 1.0: the pure biased insertion channel freezes N whenever " *
                  "the biased set is empty (every insertion becomes a null proposal and " *
                  "the matching deletions auto-reject); keep p_bias < 1 so the uniform " *
                  "sub-channel repairs ergodicity."
        end
        new(p_move, p_insert, clusters_freq, swaps_freq,
            initial_cluster_p, target_cluster_accept, cluster_adjust_interval,
            cluster_p_floor, cluster_p_ceiling,
            p_bias, bias_predicate, bias_shells, incremental)
    end
end

"""
    mutable struct GrandCanonicalNestedSamplingParameters <: SamplingParameters

Parameters for grand-canonical nested sampling on lattice systems.

The grand potential Ω = E − μN is used as the sorting quantity. Walkers have
variable particle count N, and the NS loop records (Ω, E, N) per iteration
for thermodynamic reweighting.

# Fields
- `mc_steps::Int64`: MCMC steps per replacement walker.
- `chemical_potential::Float64`: Chemical potential μ (unitless, in energy units of the Hamiltonian).
- `energy_perturbation::Float64`: Perturbation to break energy degeneracies.
- `random_seed::Int64`: Seed for the random number generator.
- `fail_count::Int64`: Consecutive failed replacements.
- `allowed_fail_count::Int64`: Maximum consecutive failures before warning.
- `init_occupation_p::Float64`: Per-site occupation probability for initial walkers.
- `n_max::Int64`: Upper bound on particle count per walker.
- `cluster_p::Float64`: Current cluster growth probability (mutable runtime state).
- `cluster_accepted::Float64`: Accepted cluster moves in current adjustment window.
- `cluster_total::Float64`: Total cluster moves attempted in current adjustment window.
- `cluster_p_history::Vector{Float64}`: Trajectory of cluster_p after each adjustment.
- `cluster_accept_history::Vector{Float64}`: Acceptance rate at each adjustment.
- `cluster_adjust_iterations::Vector{Int}`: NS iteration index at each adjustment.
- `move_stats::Dict{Symbol,Int}`: Run-total per-move-type attempt/accept counters
  accumulated from every decorrelation walk (keys match the walk's `move_stats`
  NamedTuple; cleared once at run start, never window-reset).
"""
mutable struct GrandCanonicalNestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    chemical_potential::Float64
    energy_perturbation::Float64
    random_seed::Int64
    fail_count::Int64
    allowed_fail_count::Int64
    init_occupation_p::Float64
    n_max::Int64
    cluster_p::Float64
    cluster_accepted::Float64
    cluster_total::Float64
    cluster_p_history::Vector{Float64}
    cluster_accept_history::Vector{Float64}
    cluster_adjust_iterations::Vector{Int}
    move_stats::Dict{Symbol,Int}
end

"""
    GrandCanonicalNestedSamplingParameters(;
        mc_steps=100, chemical_potential=0.0, energy_perturbation=1e-12,
        random_seed=1234, fail_count=0, allowed_fail_count=10,
        init_occupation_p=0.5, n_max=typemax(Int64),
        cluster_p=0.3, cluster_accepted=0.0, cluster_total=0.0,
        cluster_p_history=Float64[], cluster_accept_history=Float64[],
        cluster_adjust_iterations=Int[], move_stats=Dict{Symbol,Int}())

Convenience constructor for `GrandCanonicalNestedSamplingParameters`.

The `n_max` parameter sets an upper bound on the number of particles per walker.
Insertions are rejected when N ≥ n_max. Default is `typemax(Int64)` (no cap).

The `cluster_*` fields are mutable runtime state for adaptive cluster move tuning.
They are initialized from the static configuration on `MCGrandCanonicalMoves` at
the start of `grand_canonical_nested_sampling` when `clusters_freq > 0`.

`move_stats` holds run-total per-move-type attempt/accept counters accumulated
from every decorrelation walk (keys match the walk's `move_stats` NamedTuple;
never window-reset, cleared once at the start of each run).
"""
function GrandCanonicalNestedSamplingParameters(;
    mc_steps::Int64=100,
    chemical_potential::Float64=0.0,
    energy_perturbation::Float64=1e-12,
    random_seed::Int64=1234,
    fail_count::Int64=0,
    allowed_fail_count::Int64=10,
    init_occupation_p::Float64=0.5,
    n_max::Int64=typemax(Int64),
    cluster_p::Float64=0.3,
    cluster_accepted::Float64=0.0,
    cluster_total::Float64=0.0,
    cluster_p_history::Vector{Float64}=Float64[],
    cluster_accept_history::Vector{Float64}=Float64[],
    cluster_adjust_iterations::Vector{Int}=Int[],
    move_stats::Dict{Symbol,Int}=Dict{Symbol,Int}(),
)
    GrandCanonicalNestedSamplingParameters(
        mc_steps, chemical_potential, energy_perturbation,
        random_seed, fail_count, allowed_fail_count,
        init_occupation_p, n_max,
        cluster_p, cluster_accepted, cluster_total,
        cluster_p_history, cluster_accept_history, cluster_adjust_iterations,
        move_stats,
    )
end

"""
    sort_by_energy!(liveset::LJAtomWalkers)

Sorts the walkers in the liveset by their energy in descending order.

# Arguments
- `liveset::LJAtomWalkers`: The liveset of walkers to be sorted.

# Returns
- `liveset::LJAtomWalkers`: The sorted liveset.
"""
function sort_by_energy!(liveset::AbstractLiveSet)
    sort!(liveset.walkers, by = x -> x.energy, rev=true)
    # println("after sort ats[1].system_data.energy: ", ats[1].system_data.energy)
    return liveset
end

"""
    update_iter!(liveset::AtomWalkers)

Update the iteration count for each walker in the liveset.

# Arguments
- `liveset::AtomWalkers`: The set of walkers to update.

"""
function update_iter!(liveset::AbstractLiveSet)
    for at in liveset.walkers
        at.iter += 1
    end
end

const _RESERVED_LEDGER_COLUMNS = (:iter, :emax, :omega, :energy, :num_particles, :log_compression,
                                  :move_attempted, :move_accepted, :insert_attempted,
                                  :insert_accepted, :delete_attempted, :delete_accepted)

"""
    _validate_observables(observables, liveset::AbstractLiveSet)

Validate an `observables` specification — a vector of `name::Symbol =>
callback` pairs — for per-dead-point recording in a nested-sampling loop.

Throws `ArgumentError` on an empty list, duplicate names, a name colliding
with a reserved ledger column, or a callback whose probe evaluation on
`liveset.walkers[1].configuration` does not return a `Real`; warns when the
probe value is non-finite (a non-finite observable contributes NaN/Inf to
every weighted average downstream). Callbacks must be pure functions of the
configuration: each is evaluated once here as a probe and once per accepted
iteration on the culled walker.
"""
function _validate_observables(observables::AbstractVector{<:Pair{Symbol,<:Any}},
                               liveset::AbstractLiveSet)
    isempty(observables) && throw(ArgumentError(
        "observables: empty list; pass `nothing` to disable observable recording"))
    names = first.(observables)
    allunique(names) || throw(ArgumentError(
        "observables: duplicate names in $(names)"))
    for name in names
        if name in _RESERVED_LEDGER_COLUMNS
            throw(ArgumentError(
                "observables: name :$name collides with a reserved ledger " *
                "column; reserved names are $(_RESERVED_LEDGER_COLUMNS)"))
        end
    end
    for (name, f) in observables
        probe = f(liveset.walkers[1].configuration)
        probe isa Real || throw(ArgumentError(
            "observables: callback :$name returned a $(typeof(probe)); " *
            "callbacks must return a Real"))
        if !isfinite(Float64(probe))
            @warn "observables: probe evaluation of :$name returned a non-finite value ($probe)"
        end
    end
    return nothing
end

"""
    estimate_temperature(n_walker::Int, n_cull::Int, ediff::Float64)
Estimate the temperature for the nested sampling algorithm from dlog(ω)/dE.
"""
function estimate_temperature(n_walkers::Int, n_cull::Int, ediff::Float64, iter::Int=1)
    ω = (n_cull / (n_walkers + n_cull)) * (n_walkers / (n_walkers + n_cull))^iter
    β = log(ω) / ediff
    kb = 8.617333262145e-5 # eV/K
    T = 1 / (kb * β) # in Kelvin
    return T
end


"""
    _tie_block_length(walkers)

Return the number of leading walkers in an energy-sorted (descending) walker vector whose
energies are bit-exactly equal to the worst walker's energy. A return value of `1` means
the energy ceiling is unique; a larger value means the ceiling sits on an exact energy
plateau (the generic outcome of truncated pair potentials over configuration spaces with
vacuum, where the potential is exactly zero beyond every cutoff sphere).
"""
function _tie_block_length(walkers)
    e1 = walkers[1].energy
    n = 1
    while n < length(walkers) && walkers[n+1].energy == e1
        n += 1
    end
    return n
end

"""
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine)

Perform a single step of the nested sampling algorithm using the Monte Carlo random walk routine.

Exact energy ties at the ceiling (an energy plateau) are handled with plateau-aware
compression following Fowlie, Handley and Su, Mon. Not. R. Astron. Soc. 503, 1199 (2021):
when two or more live walkers tie the ceiling bit-exactly, the tied walkers are evicted
one by one without replacement, each eviction compressing the prior volume by
`(n_live - 1)/n_live` with the shrinking live count, and the live set is refilled by
cloning and decorrelating survivors below the plateau only once the last tied walker has
been evicted. A normal (unique-ceiling) cull compresses by `n_live/(n_live + 1)` as
before. The per-cull log-compression is returned as a fifth value and recorded by the
[`nested_sampling`](@ref) driver in a `log_compression` ledger column, consumed by the
log-compression method of `ωᵢ`; for tie-free ledgers the column is uniformly
`log(K/(K+1))` and the legacy iteration-based weights are unchanged. One documented
corner keeps the previous semantics: if the ENTIRE live set ties (every walker on the
plateau), no survivor samples the sub-plateau region, so the step falls back to the
ordinary clone-and-walk cull with replacement, charging `n_live/(n_live + 1)`; the
plateau under-compression bias of that corner is confined to runs whose live set is
entirely on a plateau, which for an i.i.d. initialization occurs with probability
`f^K` for plateau prior fraction `f` (about 1% at `f = 0.91`, `K = 48`, but about 21%
at `K = 16` — raise `K` when the vacuum fraction is large). This plateau handling
applies to the two SERIAL `MCRoutine` step methods (`AtomWalkers` and
`LJSurfaceWalkers`); the parallel, distributed, and mixed-moves step methods keep the
previous fixed-compression semantics and their ledgers carry no `log_compression`
column.

# Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRoutine`: The Monte Carlo routine for generating new samples. See [`MCRoutine`](@ref).

# Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
- `log_t`: The log of the prior-volume compression factor charged for this cull
  (`missing` when no walker was culled).
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    # Validate the routine before any walker is touched: the plateau branch below
    # evicts without running a walk, so an unsupported routine or dimension must
    # error here rather than silently evicting on tie steps.
    mc_routine isa Union{MCRandomWalkMaxE, MCRandomWalkClone} || error("Unsupported MCRoutine type: $mc_routine")
    length(mc_routine.dims) in (2, 3) || error("Unsupported dimensions: $(mc_routine.dims)")
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    log_t::Union{Missing,Float64} = missing
    n_live = length(ats)
    n_tied = _tie_block_length(ats)
    if (n_tied >= 2 || ns_params.plateau_refill_target != 0) && n_tied < n_live
        # Plateau block: evict the worst tied walker without replacement
        # (Fowlie-Handley-Su), charging (n_live - 1)/n_live with the shrinking
        # live count. With the refill target set and a unique ceiling, this is
        # the last tied walker of the block: evict it, then refill.
        if ns_params.plateau_refill_target == 0
            ns_params.plateau_refill_target = n_live
        end
        popfirst!(ats)
        update_iter!(liveset)
        iter = liveset.walkers[1].iter
        log_t = log((n_live - 1) / n_live)
        if n_tied == 1
            # The plateau is exhausted: refill by cloning survivors and
            # decorrelating strictly below the plateau energy (volume-neutral;
            # no ledger row), then clear the block state.
            refill_fails = 0
            while length(ats) < ns_params.plateau_refill_target && refill_fails < ns_params.allowed_fail_count
                # The refill walk uses the routine's own dimension constraint so
                # refilled clones stay on the constrained prior manifold.
                if length(mc_routine.dims) == 3
                    accept_r, _, at_r = MC_random_walk!(ns_params.mc_steps, deepcopy(rand(ats)), lj, ns_params.step_size, emax)
                else
                    accept_r, _, at_r = MC_random_walk_2D!(ns_params.mc_steps, deepcopy(rand(ats)), lj, ns_params.step_size, emax; dims=mc_routine.dims)
                end
                if accept_r
                    push!(ats, at_r)
                else
                    refill_fails += 1
                end
            end
            length(ats) < ns_params.plateau_refill_target && @warn "Plateau refill left the live set at $(length(ats)) of $(ns_params.plateau_refill_target) walkers after $(refill_fails) failed decorrelation attempts; subsequent culls are charged with the actual live count."
            ns_params.plateau_refill_target = 0
        end
        return iter, emax, liveset, ns_params, log_t
    end
    if mc_routine isa MCRandomWalkMaxE
        to_walk = deepcopy(ats[1])
    elseif mc_routine isa MCRandomWalkClone
        to_walk = deepcopy(rand(ats[2:end]))
    else
        error("Unsupported MCRoutine type: $mc_routine")
    end
    if length(mc_routine.dims) == 3
        accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    elseif length(mc_routine.dims) == 2
        accept, rate, at = MC_random_walk_2D!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax; dims=mc_routine.dims)
        # @info "Doing a 2D random walk"
    elseif length(mc_routine.dims) == 1
        error("Unsupported dimensions: $(mc_routine.dims)")
    end
    # accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: $(round(ns_params.step_size; sigdigits=4))"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
        log_t = log(n_live / (n_live + 1))
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax, liveset, ns_params, log_t
end

"""
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutineParallel)
Perform a single step of the nested sampling algorithm using the parallel Monte Carlo random walk routine.
# Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRoutineParallel`: The parallel Monte Carlo routine for generating new samples. See [`MCRoutineParallel`](@ref).
# Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDistributed; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    to_walk_inds = sample(2:length(ats), nworkers(); replace=false)
    # println("to_walk_inds: ", to_walk_inds) # DEBUG
    
    to_walks = deepcopy.(ats[to_walk_inds])

    if length(mc_routine.dims) == 3
        random_walk_function = MC_random_walk!
        walk_kwargs = NamedTuple()
    elseif length(mc_routine.dims) == 2
        random_walk_function = MC_random_walk_2D!
        walk_kwargs = (dims = mc_routine.dims,)
    else
        error("Unsupported dimensions: $(mc_routine.dims)")
    end

    mc_steps_per_worker = ceil(Int, ns_params.mc_steps / nworkers()) # distribute the total MC steps among workers

    walking = [remotecall(random_walk_function, workers()[i], mc_steps_per_worker, to_walk, lj, ns_params.step_size, emax[mc_routine.n_cull]; walk_kwargs...) for (i,to_walk) in enumerate(to_walks)]
    walked = fetch.(walking)
    finalize.(walking) # finalize the remote calls, clear the memory

    accepted_rates = [x[2] for x in walked]
    rate = mean(accepted_rates)

    # sort!(walked, by = x -> x[3].energy, rev=true)
    # filter!(x -> x[1], walked) # remove the failed ones
    accepted_inds = findall(x -> x[1]==1, walked)

    if length(accepted_inds) < mc_routine.n_cull # if not enough accepted walkers
        ns_params.fail_count += 1
        emax = [missing]
        return iter, emax[end], liveset, ns_params
    else
        # pick one from the accepted ones
        picked = sample(accepted_inds, mc_routine.n_cull; replace=false)
        for (i, ind) in enumerate(picked)
            ats[i] = walked[ind][3]
        end
        # println("picked: ", picked) # DEBUG
        # remove the picked one from accepted_inds
        filter!(x -> x ∉ picked, accepted_inds)
        # println("remaining accepted_inds: ", accepted_inds) # DEBUG

        if !isempty(accepted_inds)
            for i in accepted_inds
                ats[to_walk_inds[i]] = walked[i][3]
                # println("Updating ats at index $(to_walk_inds[i])") # DEBUG
            end
        end
    end

    update_iter!(liveset)
    ns_params.fail_count = 0
    iter = liveset.walkers[1].iter

    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax[mc_routine.n_cull], liveset, ns_params
end

function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutineParallel; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    if mc_routine isa MCRandomWalkMaxEParallel
        to_walk_inds = 1:nworkers()
    elseif mc_routine isa MCRandomWalkCloneParallel
        to_walk_inds = sample(2:length(ats), nworkers(); replace=false)
    end
    
    to_walks = deepcopy.(ats[to_walk_inds])

    if length(mc_routine.dims) == 3
        random_walk_function = MC_random_walk!
        walk_kwargs = NamedTuple()
    elseif length(mc_routine.dims) == 2
        random_walk_function = MC_random_walk_2D!
        walk_kwargs = (dims = mc_routine.dims,)
    else
        error("Unsupported dimensions: $(mc_routine.dims)")
    end


    walking = [remotecall(random_walk_function, workers()[i], ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax[end]; walk_kwargs...) for (i,to_walk) in enumerate(to_walks)]
    walked = fetch.(walking)
    finalize.(walking) # finalize the remote calls, clear the memory

    accepted_rates = [x[2] for x in walked]
    rate = mean(accepted_rates)

    if prod([x[1] for x in walked]) == 0 # if any of the walkers failed
        ns_params.fail_count += 1
        emax = [missing]
        return iter, emax[end], liveset, ns_params
    end

    # sort!(walked, by = x -> x[3].energy, rev=true)
    # filter!(x -> x[1], walked) # remove the failed ones

    for (i, at) in enumerate(walked)
        ats[i] = at[3]
    end

    update_iter!(liveset)
    ns_params.fail_count = 0
    iter = liveset.walkers[1].iter

    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax[end], liveset, ns_params
end

function nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutineParallel; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    if mc_routine isa MCRandomWalkMaxEParallel
        to_walk_inds = 1:nworkers()
    elseif mc_routine isa MCRandomWalkCloneParallel
        to_walk_inds = sample(2:length(ats), nworkers(); replace=false)
    end
    
    to_walks = deepcopy.(ats[to_walk_inds])

    if length(mc_routine.dims) == 3
        random_walk_function = MC_random_walk!
    else
        # The surface walks are 3D-only, matching the serial surface step's
        # restriction: a length-2 dims would select a walk method with no
        # surface variant and fail as a MethodError on the worker.
        error("Unsupported dimensions: $(mc_routine.dims)")
    end


    walking = [remotecall(random_walk_function, workers()[i], ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax[end], liveset.surface) for (i,to_walk) in enumerate(to_walks)]
    walked = fetch.(walking)
    finalize.(walking) # finalize the remote calls, clear the memory

    accepted_rates = [x[2] for x in walked]
    rate = mean(accepted_rates)

    if prod([x[1] for x in walked]) == 0 # if any of the walkers failed
        ns_params.fail_count += 1
        emax = [missing]
        return iter, emax[end], liveset, ns_params
    end

    # sort!(walked, by = x -> x[3].energy, rev=true)
    # filter!(x -> x[1], walked) # remove the failed ones

    for (i, at) in enumerate(walked)
        ats[i] = at[3]
    end

    update_iter!(liveset)
    ns_params.fail_count = 0
    iter = liveset.walkers[1].iter

    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax[end], liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine)

Serial nested-sampling step for surface livesets. Identical to the `AtomWalkers` serial
method — including the plateau-aware handling of exact energy ties, the fifth
`log_t` return value, and the all-tied fallback documented there — except that every
decorrelation and refill walk carries the frozen surface and only `dims == [1, 2, 3]`
is supported.
"""
function nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    # Validate the routine before any walker is touched (the plateau branch below
    # evicts without running a walk); this method's walks are 3D-only.
    mc_routine isa Union{MCRandomWalkMaxE, MCRandomWalkClone} || error("Unsupported MCRoutine type: $mc_routine")
    length(mc_routine.dims) == 3 || error("Unsupported dimensions: $(mc_routine.dims)")
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    log_t::Union{Missing,Float64} = missing
    n_live = length(ats)
    n_tied = _tie_block_length(ats)
    if (n_tied >= 2 || ns_params.plateau_refill_target != 0) && n_tied < n_live
        # Plateau block: evict the worst tied walker without replacement
        # (Fowlie-Handley-Su); see the AtomWalkers method docstring.
        if ns_params.plateau_refill_target == 0
            ns_params.plateau_refill_target = n_live
        end
        popfirst!(ats)
        update_iter!(liveset)
        iter = liveset.walkers[1].iter
        log_t = log((n_live - 1) / n_live)
        if n_tied == 1
            refill_fails = 0
            while length(ats) < ns_params.plateau_refill_target && refill_fails < ns_params.allowed_fail_count
                accept_r, _, at_r = MC_random_walk!(ns_params.mc_steps, deepcopy(rand(ats)), lj, ns_params.step_size, emax, liveset.surface)
                if accept_r
                    push!(ats, at_r)
                else
                    refill_fails += 1
                end
            end
            length(ats) < ns_params.plateau_refill_target && @warn "Plateau refill left the live set at $(length(ats)) of $(ns_params.plateau_refill_target) walkers after $(refill_fails) failed decorrelation attempts; subsequent culls are charged with the actual live count."
            ns_params.plateau_refill_target = 0
        end
        return iter, emax, liveset, ns_params, log_t
    end
    if mc_routine isa MCRandomWalkMaxE
        to_walk = deepcopy(ats[1])
    elseif mc_routine isa MCRandomWalkClone
        to_walk = deepcopy(rand(ats[2:end]))
    else
        error("Unsupported MCRoutine type: $mc_routine")
    end
    if length(mc_routine.dims) == 3
        accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax, liveset.surface)
    else
        error("Unsupported dimensions: $(mc_routine.dims)")
    end
    # accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: $(round(ns_params.step_size; sigdigits=4))"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
        log_t = log(n_live / (n_live + 1))
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax, liveset, ns_params, log_t
end

"""
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCMixedMoves)

Perform a single step of the nested sampling algorithm using the Monte Carlo mixed moves routine.
By default, this routine performs parallel decorrelation of multiple walkers.

Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCMixedMoves`: The Monte Carlo mixed moves routine.

Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.

Note
- To invoke the parallel version of this routine, use `MCMixedMovesParallel` as the `mc_routine` argument.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCMixedMoves; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    to_walk = deepcopy(rand(ats[2:end]))

    accept, rate, at = MC_mixed_moves!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax, [mc_routine.walks_freq, mc_routine.swaps_freq])

    # accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: $(round(ns_params.step_size; sigdigits=4))"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate; range=ns_params.accept_range)

    return iter, emax, liveset, ns_params
end

function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCMixedMovesParallel; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    to_walk_inds = sample(2:length(ats), nworkers(); replace=false)
    # println("to_walk_inds: ", to_walk_inds) # DEBUG
    
    to_walks = deepcopy.(ats[to_walk_inds])

    walking = [remotecall(MC_mixed_moves!, workers()[i], ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax[1], [mc_routine.walks_freq, mc_routine.swaps_freq]) for (i,to_walk) in enumerate(to_walks)]
    walked = fetch.(walking)
    finalize.(walking) # finalize the remote calls, clear the memory

    accepted_rates = [x[2] for x in walked]
    rate = mean(accepted_rates)

    # sort!(walked, by = x -> x[3].energy, rev=true)
    # filter!(x -> x[1], walked) # remove the failed ones
    accepted_inds = findall(x -> x[1]==1, walked)

    if length(accepted_inds) == 0 # if all of the walkers failed
        ns_params.fail_count += 1
        emax = [missing]
        return iter, emax[end], liveset, ns_params
    else
        # pick one from the accepted ones
        picked = rand(accepted_inds)
        ats[1] = walked[picked][3]
        # println("picked: ", picked) # DEBUG
        # remove the picked one from accepted_inds
        filter!(x -> x != picked, accepted_inds)
        # println("remaining accepted_inds: ", accepted_inds) # DEBUG

        if !isempty(accepted_inds)
            for i in accepted_inds
                ats[to_walk_inds[i]] = walked[i][3]
                # println("Updating ats at index $(to_walk_inds[i])") # DEBUG
            end
        end
    end

    update_iter!(liveset)
    ns_params.fail_count = 0
    iter = liveset.walkers[1].iter

    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax[1], liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers, ns_params::NestedSamplingParameters, mc_routine::MCMixedMoves)

Perform a single step of the nested sampling algorithm using a mix of geometric cluster moves and local swap moves.

The total `mc_steps` from `ns_params` are split between cluster moves and local swaps according to `clusters_freq` and
`walks_freq` in the `mc_routine`. Cluster moves use `geometric_cluster_swap!` with growth probability `ns_params.cluster_p`,
which is adaptively tuned to maintain `mc_routine.target_cluster_accept`. Local swaps use the standard `lattice_random_walk!`.

## Arguments
- `liveset::LatticeGasWalkers`: The liveset of lattice gas walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCMixedMoves`: The mixed moves routine with cluster and local swap frequencies.

## Returns
- `iter`: The iteration number of the liveset after the step.
- `emax`: The maximum energy of the liveset after the step.
- `liveset::LatticeGasWalkers`: The updated liveset.
- `ns_params::NestedSamplingParameters`: The updated parameters.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               ns_params::NestedSamplingParameters,
                               mc_routine::MCMixedMoves;
                               ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,Float64} = liveset.walkers[1].energy.val

    # Clone a random non-worst walker
    to_walk = deepcopy(rand(ats[2:end]))

    # Compute move counts from frequencies
    total_freq = mc_routine.walks_freq + mc_routine.clusters_freq
    n_local = round(Int, ns_params.mc_steps * mc_routine.walks_freq / max(total_freq, 1))
    n_cluster = ns_params.mc_steps - n_local

    # Apply cluster moves
    cluster_accepted = false
    cluster_rate = 0.0
    if n_cluster > 0
        cluster_accepted, cluster_rate, to_walk = MC_cluster_walk!(
            n_cluster, to_walk, h, emax, ns_params.cluster_p;
            energy_perturb=ns_params.energy_perturbation)
    end

    # Apply local swap moves
    local_accepted = false
    local_rate = 0.0
    if n_local > 0
        local_accepted, local_rate, to_walk = MC_random_walk!(
            n_local, to_walk, h, emax;
            energy_perturb=ns_params.energy_perturbation)
    end

    accept = cluster_accepted || local_accepted
    if accept
        push!(ats, to_walk)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        emax = missing
        ns_params.fail_count += 1
    end

    # Accumulate cluster acceptance stats for adaptive tuning
    if n_cluster > 0
        ns_params.cluster_accepted += cluster_rate * n_cluster
        ns_params.cluster_total += n_cluster
        if mc_routine.cluster_adjust_interval > 0 &&
           ns_params.cluster_total >= mc_routine.cluster_adjust_interval * n_cluster
            window_rate = ns_params.cluster_accepted / max(ns_params.cluster_total, 1.0)
            adjust_cluster_p(ns_params, window_rate, ns_iteration;
                             target=mc_routine.target_cluster_accept,
                             floor=mc_routine.cluster_p_floor,
                             ceiling=mc_routine.cluster_p_ceiling)
            ns_params.cluster_accepted = 0.0
            ns_params.cluster_total = 0.0
        end
    end

    return iter, emax * unit(liveset.walkers[1].energy), liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers, ns_params::LatticeNestedSamplingParameters, mc_routine::MCRoutine)

Perform a single step of the nested sampling algorithm.

This function takes a `liveset` of lattice gas walkers, `ns_params` containing the parameters for nested sampling, and `mc_routine` representing the Monte Carlo
routine for generating new samples. It performs a single step of the nested sampling algorithm by updating the liveset with a new walker.

## Arguments
- `liveset::LatticeGasWalkers`: The liveset of lattice gas walkers.
- `ns_params::LatticeNestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRoutine`: The Monte Carlo routine for generating new samples.

## Returns
- `iter`: The iteration number of the liveset after the step.
- `emax`: The maximum energy of the liveset after the step.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               ns_params::NestedSamplingParameters,
                               mc_routine::MCRoutine;
                               ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,Float64} = liveset.walkers[1].energy.val
    if mc_routine isa MCRandomWalkMaxE
        to_walk = deepcopy(ats[1])
    elseif mc_routine isa MCRandomWalkClone
        to_walk = deepcopy(rand(ats[2:end]))
    else
        error("Unsupported MCRoutine type: $mc_routine")
    end
    accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, h, emax; energy_perturb=ns_params.energy_perturbation)

    # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $rate, emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    # adjust_step_size(ns_params, rate)
    return iter, emax * unit(liveset.walkers[1].energy), liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers, ns_params::LatticeNestedSamplingParameters, mc_routine::MCNewSample)

Perform a single step of the nested sampling algorithm.

This function takes a `liveset` of lattice gas walkers, `ns_params` containing the parameters for nested sampling, and `mc_routine` representing the Monte Carlo routine for generating new samples. It performs a single step of the nested sampling algorithm by updating the liveset with a new walker.

## Arguments
- `liveset::LatticeGasWalkers`: The liveset of lattice gas walkers.
- `ns_params::LatticeNestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCNewSample`: The Monte Carlo routine for generating new samples.

## Returns
- `iter`: The iteration number of the liveset after the step.
- `emax`: The maximum energy of the liveset after the step.
- `liveset::LatticeGasWalkers`: The updated liveset after the step.
- `ns_params::LatticeNestedSamplingParameters`: The updated nested sampling parameters after the step.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               ns_params::NestedSamplingParameters,
                               mc_routine::MCNewSample;
                               ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,Float64} = liveset.walkers[1].energy.val

    to_walk = deepcopy(ats[1])

    accept, at = MC_new_sample!(to_walk, h, emax; energy_perturb=ns_params.energy_perturbation)

    # @info "iter: $(liveset.walkers[1].iter), emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    # adjust_step_size(ns_params, rate)
    return iter, emax * unit(liveset.walkers[1].energy), liveset, ns_params
end


function nested_sampling_step!(liveset::LatticeGasWalkers,
                               ns_params::NestedSamplingParameters,
                               mc_routine::MCRejectionSampling;
                               ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,Float64} = liveset.walkers[1].energy.val

    to_walk = deepcopy(ats[1])

    accept, at = MC_rejection_sampling!(to_walk, h, emax; energy_perturb=ns_params.energy_perturbation)

    # @info "iter: $(liveset.walkers[1].iter), emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    # adjust_step_size(ns_params, rate)
    return iter, emax * unit(liveset.walkers[1].energy), liveset, ns_params
end



"""
    nested_sampling(liveset::AbstractLiveSet, ns_params::NestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine; args...)

Perform a nested sampling loop for a given number of steps.

# Arguments
- `liveset::AbstractLiveSet`: The initial set of walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `n_steps::Int64`: The number of steps to perform.
- `mc_routine::MCRoutine`: The Monte Carlo routine to use.
- `observables`: Optional vector of `name::Symbol => callback` pairs. Each
  callback is evaluated on the culled walker's `configuration` at every
  accepted iteration and recorded as an extra `Float64` ledger column named
  `name`, exactly paired with that row's `iter`-keyed prior-volume weight.
  Names must not collide with the built-in ledger columns
  (`$(join(String.(_RESERVED_LEDGER_COLUMNS), ", "))`); callbacks must be
  pure functions of the configuration returning a `Real`. Parallel and
  multi-cull MC routines (`MCRoutineParallel` subtypes) are rejected up
  front with an `ArgumentError`; a bit-exact pairing guard additionally
  raises an error if a step ever culls a different walker than the one the
  row records. The default `nothing` leaves the ledger schema unchanged.
- `dead_point_callback`: Optional function called once per recorded dead
  point as `dead_point_callback(iter, walker)`, immediately after the
  ledger row is pushed, with the same culled walker the row describes
  (protected by the same bit-exact pairing guard as `observables`, which is
  engaged whenever either keyword is supplied). Iterations that record no
  dead point invoke no callback, so the invocation count equals `nrow(df)`.
  The walker is the loop's live object: treat it as read-only and copy
  anything kept, e.g. `copy(walker.configuration.components[1])` or
  `deepcopy(walker.configuration)`. Composes freely with `observables`; the
  surviving live set needs no callback because the loop returns it.
  Parallel and multi-cull MC routines are rejected up front. Typical uses:
  collecting occupation vectors for batch post-processing, and re-evaluating
  each dead point under a second Hamiltonian, whose energies substitute for
  the ledger energy column in the grand-canonical stats functions (the
  prior-volume weights depend only on the iteration index, so the
  substituted column is a valid estimator over the same ladder). The
  default `nothing` changes nothing.

# Returns
- `df`: A DataFrame containing the iteration number and maximum energy for each step,
  plus one column per requested observable.
- `liveset`: The updated set of walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling(liveset::AbstractLiveSet,
                                ns_params::NestedSamplingParameters,
                                n_steps::Int64,
                                mc_routine::MCRoutine,
                                save_strategy::DataSavingStrategy;
                                observables::Union{Nothing,AbstractVector{<:Pair{Symbol,<:Any}}}=nothing,
                                dead_point_callback::Union{Nothing,Function}=nothing)
    # Initialize cluster_p and reset counters from MCMixedMoves if applicable
    if mc_routine isa MCMixedMoves && mc_routine.clusters_freq > 0
        ns_params.cluster_p = mc_routine.initial_cluster_p
        ns_params.cluster_accepted = 0.0
        ns_params.cluster_total = 0.0
        empty!(ns_params.cluster_p_history)
        empty!(ns_params.cluster_accept_history)
        empty!(ns_params.cluster_adjust_iterations)
    end
    # Defensive: a parameters object recovered from a run that ended inside a
    # plateau block must not arm a fresh run's first cull as a tie eviction.
    ns_params.plateau_refill_target = 0
    df = DataFrame(iter=Int[], emax=Float64[])
    if observables !== nothing
        mc_routine isa MCRoutineParallel && throw(ArgumentError(
            "observables: parallel/multi-cull MC routines are not supported " *
            "— each accepted iteration culls multiple walkers but records a " *
            "single ledger row, so per-dead-point pairing is undefined"))
        _validate_observables(observables, liveset)
        for (name, _) in observables
            df[!, name] = Float64[]
        end
    end
    if dead_point_callback !== nothing
        mc_routine isa MCRoutineParallel && throw(ArgumentError(
            "dead_point_callback: parallel/multi-cull MC routines are not " *
            "supported — each accepted iteration culls multiple walkers but " *
            "invokes a single callback, so per-dead-point pairing is " *
            "undefined"))
    end
    culled = nothing
    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)
        if observables !== nothing || dead_point_callback !== nothing
            # Pre-sort with the step's own comparator and hold the walker the
            # step will cull: the step re-sorts (stable sort, so the identity
            # permutation here) and culls walkers[1] on accept. `popfirst!`
            # removes but never mutates the culled walker, so evaluating the
            # callbacks on it after the step is exact.
            sort_by_energy!(liveset)
            culled = liveset.walkers[1]
        end
        step_ret = nested_sampling_step!(liveset, ns_params, mc_routine; ns_iteration=i)
        iter, emax, liveset, ns_params = step_ret[1], step_ret[2], step_ret[3], step_ret[4]
        # Serial atomistic steps additionally return the per-cull log-compression
        # (plateau-aware culling); other step methods keep the four-value return
        # and their ledgers gain no column.
        log_t = length(step_ret) >= 5 ? step_ret[5] : missing
        @debug "n_step $i, iter: $iter, emax: $emax"
        if ns_params.fail_count >= ns_params.allowed_fail_count
            @warn "Failed to accept MC move $(ns_params.allowed_fail_count) times in a row. Reset step size!"
            ns_params.fail_count = 0
            ns_params.step_size = ns_params.initial_step_size
        end
        if !(iter isa typeof(missing))
            if culled !== nothing
                emax == culled.energy || error(
                    "NS observable recording: dead-point/observable pairing " *
                    "lost (step culled a walker with energy $emax, but the " *
                    "pre-sorted worst walker had $(culled.energy)); this MC " *
                    "routine is not supported with `observables`/" *
                    "`dead_point_callback`")
            end
            if log_t !== missing && !hasproperty(df, :log_compression)
                # Every accepted row of a log-compression-emitting run carries a
                # value, and the first accepted row is the first push, so the
                # column can be added lazily without backfilling.
                df[!, :log_compression] = Float64[]
            end
            row = observables === nothing ? (iter, emax.val) :
                (iter, emax.val,
                 (Float64(f(culled.configuration)) for (_, f) in observables)...)
            push!(df, log_t === missing ? row : (row..., log_t))
            dead_point_callback === nothing || dead_point_callback(iter, culled)
        end
        print_message(i, iter, emax, ns_params.step_size, print_info, liveset)
        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end
    return df, liveset, ns_params
end

function print_message(i, iter, emax, step_size, print_info, liveset::LatticeWalkers)
    if print_info && !(iter isa typeof(missing))
        @info "iter: $(liveset.walkers[1].iter), emax: $(emax)"
    elseif print_info && iter isa typeof(missing)
        @info "MC move failed, step: $(i), emax: $(liveset.walkers[1].energy)"
    end
end

function print_message(i, iter, emax, step_size, print_info, liveset::AtomWalkers)
    if print_info && !(iter isa typeof(missing))
        @info "iter: $(liveset.walkers[1].iter), emax: $(emax.val), step_size: $(round(step_size; sigdigits=4))"
    elseif print_info && iter isa typeof(missing)
        @info "MC move failed, step: $(i), emax: $(liveset.walkers[1].energy.val), step_size: $(round(step_size; sigdigits=4))"
    end
end


# ======================================================================
# Grand-canonical nested sampling
# ======================================================================

"""
    _grand_potential(walker::LatticeWalker{1}, mu::Float64) -> typeof(0.0u"eV")

Compute Ω = E − μN for a single-component lattice walker.
"""
function _grand_potential(walker::LatticeWalker{1}, mu::Float64)
    n = sum(walker.configuration.components[1])
    return walker.energy - mu * n * unit(walker.energy)
end

"""
    _clone_walker_shared_geometry(w::LatticeWalker)

Clone a lattice walker for a replacement walk, copying only its occupancy
vectors, energy, and iteration counter while sharing the run-invariant
geometry by reference (the `Val(:share_geometry)` constructor). The
geometry fields are never mutated after construction and the lattice
drivers are serial loops, so the clone is behaviorally identical to a
`deepcopy`: it draws no randomness and changes no floating-point value.
"""
_clone_walker_shared_geometry(w::LatticeWalker) =
    _shared_geometry_walker(w, w.configuration)

function _shared_geometry_walker(w::LatticeWalker, cfg::MLattice{C,G}) where {C,G}
    shared = MLattice{C,G}(Val(:share_geometry), cfg,
                           [copy(v) for v in cfg.components])
    return LatticeWalker(shared, energy=w.energy, iter=w.iter)
end

"""
    _perturbation_energy_bound(h, lattice) -> Union{Float64,Nothing}

Magnitude bound on the lattice energy in the Hamiltonian's own energy
units: the on-site magnitude times the site count, plus each coupled
shell's magnitude times its total ordered-entry count over two (image
multiplicity included, read from the built neighbor lists), plus the
site-field and cluster magnitudes where applicable. Returns `nothing` for
Hamiltonian types without a method, which skips the tie-breaker warning.
"""
function _perturbation_energy_bound(h::GenericLatticeHamiltonian{N,U},
                                    lattice) where {N,U}
    nbrs = lattice.neighbors
    b = abs(ustrip(h.on_site_interaction)) * length(lattice.components[1])
    shells = isempty(nbrs) ? 0 : length(nbrs[1])
    for n in 1:min(N, shells)
        entries = sum(length(nbrs[i][n]) for i in eachindex(nbrs))
        b += abs(ustrip(h.nth_neighbor_interactions[n])) * entries / 2
    end
    return b
end
function _perturbation_energy_bound(h::SiteFieldLatticeHamiltonian, lattice)
    b = _perturbation_energy_bound(h.base, lattice)
    b === nothing && return nothing
    return b + sum(x -> abs(ustrip(x)), h.field)
end
_perturbation_energy_bound(h::MLatticeHamiltonian, lattice) =
    _perturbation_energy_bound(h.Hamiltonians[1, 1], lattice)
function _perturbation_energy_bound(h::ClusterLatticeHamiltonian, lattice)
    b = _perturbation_energy_bound(h.pair_ham, lattice)
    b === nothing && return nothing
    for c in h.clusters
        b += abs(ustrip(c.coupling)) * length(c.embeddings)
    end
    return b
end
_perturbation_energy_bound(h, lattice) = nothing

"""
    _warn_perturbation_scale(liveset::LatticeGasWalkers, delta::Float64)

Warn once, at driver entry, when the configured `energy_perturbation` sits
below the documented degenerate-plateau resolution floor
`K^2 * eps(E_bound)`: below it the K walkers cannot reliably draw
distinct tie-breaking values on a degenerate plateau, and the failure mode
is silent (plateaus decompress at machine granularity and read as poor
mixing). A warning, never a throw; a zero perturbation is a deliberate
no-perturbation choice and is not warned about.
"""
function _warn_perturbation_scale(liveset::LatticeGasWalkers, delta::Float64)
    isempty(liveset.walkers) && return nothing
    bound = _perturbation_energy_bound(liveset.hamiltonian,
                                       liveset.walkers[1].configuration)
    bound === nothing && return nothing
    floor_delta = length(liveset.walkers)^2 * eps(bound)
    if 0.0 < delta < floor_delta
        @warn "energy_perturbation = $delta sits below the degenerate-" *
              "plateau resolution floor K^2 * eps(E_bound) = $floor_delta " *
              "(E_bound = $bound in the Hamiltonian's energy units, K = " *
              "$(length(liveset.walkers))); ties may decompress at machine " *
              "granularity and read as poor mixing. Recommended minimum: " *
              "$floor_delta."
    end
    return nothing
end

"""
    _init_gc_walkers!(liveset::LatticeGasWalkers, gc_params::GrandCanonicalNestedSamplingParameters)

Initialize walkers with random microstates for grand-canonical NS.
Each site is occupied independently with probability `gc_params.init_occupation_p`.
"""
function _init_gc_walkers!(liveset::LatticeGasWalkers, gc_params::GrandCanonicalNestedSamplingParameters)
    h = liveset.hamiltonian
    n_max = gc_params.n_max
    for walker in liveset.walkers
        random_microstate!(walker.configuration; p=gc_params.init_occupation_p)
        # Enforce n_max: if too many particles, randomly delete until N ≤ n_max
        n_occ = sum(walker.configuration.components[1])
        if n_occ > n_max
            occupied = findall(walker.configuration.components[1])
            shuffle!(occupied)
            for i in 1:(n_occ - n_max)
                walker.configuration.components[1][occupied[i]] = false
            end
        end
        assign_energy!(walker, h; perturb_energy=gc_params.energy_perturbation)
    end
    return liveset
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers,
                          gc_params::GrandCanonicalNestedSamplingParameters,
                          mc_routine::MCGrandCanonicalMoves;
                          ns_iteration::Int=0)

Perform one step of grand-canonical nested sampling.

Sorts walkers by Ω = E − μN, removes the worst (highest Ω), clones a parent
with Ω < Ω_worst, and decorrelates the clone via grand-canonical MCMC.

# Returns
- `iter`: Iteration number (or `missing` if the step failed).
- `omega_max`: The Ω value of the removed walker (with units).
- `energy`: The E value of the removed walker.
- `num_particles`: The N value of the removed walker.
- `liveset`: The updated liveset.
- `gc_params`: The updated parameters.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               gc_params::GrandCanonicalNestedSamplingParameters,
                               mc_routine::MCGrandCanonicalMoves;
                               ns_iteration::Int=0)
    ats = liveset.walkers
    h = liveset.hamiltonian
    mu = gc_params.chemical_potential
    n_walkers = length(ats)

    # Sort by grand potential (descending — worst first): one O(K·M) key
    # sweep, then a stable sortperm on the cached keys. The comparisons see
    # the identical key values, so the resulting order (ties included), the
    # random stream, and every recorded digit match the by-comparator sort,
    # while the per-iteration cost drops from O(K·M·log K) to
    # O(K·M + K·log K); the cached keys are reused below.
    omega_keys = [_grand_potential(w, mu) for w in ats]
    omega_perm = sortperm(omega_keys, rev=true)
    permute!(ats, omega_perm)
    permute!(omega_keys, omega_perm)

    iter::Union{Missing,Int} = missing
    worst = ats[1]
    omega_worst = omega_keys[1]
    energy_worst = worst.energy
    n_worst = sum(worst.configuration.components[1])

    # Select parent: prefer walkers strictly below omega_worst
    omega_max_val = omega_worst.val  # unitless for the MC function
    eligible = [k for k in 2:n_walkers if omega_keys[k] < omega_worst]
    if !isempty(eligible)
        parent_idx = rand(eligible)
    else
        parent_idx = rand(2:n_walkers)
    end
    to_walk = _clone_walker_shared_geometry(ats[parent_idx])

    # Decorrelate via GC MCMC
    accept, rate, to_walk, cl_accepted, cl_total, move_stats = MC_grand_canonical_walk!(
        gc_params.mc_steps, to_walk, h, omega_max_val, mu;
        p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
        energy_perturb=gc_params.energy_perturbation,
        n_max=gc_params.n_max,
        clusters_freq=mc_routine.clusters_freq,
        swaps_freq=mc_routine.swaps_freq,
        cluster_p=gc_params.cluster_p,
        p_bias=mc_routine.p_bias,
        bias_predicate=mc_routine.bias_predicate,
        bias_shells=mc_routine.bias_shells,
        incremental=mc_routine.incremental)

    if accept
        push!(ats, to_walk)
        popfirst!(ats)
        update_iter!(liveset)
        gc_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        omega_worst = missing
        energy_worst = missing
        n_worst = missing
        gc_params.fail_count += 1
    end

    # Accumulate cluster acceptance stats and adapt cluster_p
    _accumulate_cluster_stats!(gc_params, mc_routine, cl_accepted, cl_total, ns_iteration)
    _accumulate_move_stats!(gc_params, move_stats)

    return iter, omega_worst, energy_worst, n_worst, liveset, gc_params
end

"""
    grand_canonical_nested_sampling(liveset::LatticeGasWalkers,
                                    gc_params::GrandCanonicalNestedSamplingParameters,
                                    n_steps::Int64,
                                    mc_routine::MCGrandCanonicalMoves,
                                    save_strategy::DataSavingStrategy)

Run the grand-canonical nested sampling loop.

Initializes walkers with random microstates, then iterates: remove the
highest-Ω walker, record (Ω, E, N), replace with a decorrelated clone.

# Arguments
- `liveset::LatticeGasWalkers`: The initial liveset (walkers will be re-initialized).
- `gc_params::GrandCanonicalNestedSamplingParameters`: GC-NS parameters including μ.
- `n_steps::Int64`: Number of NS iterations.
- `mc_routine::MCGrandCanonicalMoves`: The GC move routine.
- `save_strategy::DataSavingStrategy`: Strategy for periodic output.
- `observables`: Optional vector of `name::Symbol => callback` pairs
  recording per-dead-point observables as extra ledger columns; see
  [`nested_sampling`](@ref).
- `dead_point_callback`: Optional per-dead-point callback
  `(iter, walker) -> ...` invoked with the culled walker immediately after
  its ledger row is pushed; see [`nested_sampling`](@ref) for the contract.

- `stop_on_stall::Bool=false`: When true and `fail_count` reaches
  `allowed_fail_count`, warn once and return the partial ledger and the
  intact live set (`fail_count` stays at threshold); the default keeps the
  shipped warn-and-continue behavior byte-identically.

# Returns
- `df::DataFrame`: Columns `[:iter, :omega, :energy, :num_particles]`,
  plus one column per requested observable.
- `liveset::LatticeGasWalkers`: The final liveset (surviving walkers).
- `gc_params::GrandCanonicalNestedSamplingParameters`: Updated parameters.
"""
function grand_canonical_nested_sampling(liveset::LatticeGasWalkers,
                                         gc_params::GrandCanonicalNestedSamplingParameters,
                                         n_steps::Int64,
                                         mc_routine::MCGrandCanonicalMoves,
                                         save_strategy::DataSavingStrategy;
                                         observables::Union{Nothing,AbstractVector{<:Pair{Symbol,<:Any}}}=nothing,
                                         dead_point_callback::Union{Nothing,Function}=nothing,
                                         stop_on_stall::Bool=false)
    # Initialize walkers with random microstates
    _init_gc_walkers!(liveset, gc_params)
    _warn_perturbation_scale(liveset, gc_params.energy_perturbation)

    # Initialize cluster_p and reset counters from MCGrandCanonicalMoves if applicable
    if mc_routine.clusters_freq > 0
        gc_params.cluster_p = mc_routine.initial_cluster_p
        gc_params.cluster_accepted = 0.0
        gc_params.cluster_total = 0.0
        empty!(gc_params.cluster_p_history)
        empty!(gc_params.cluster_accept_history)
        empty!(gc_params.cluster_adjust_iterations)
    end
    # Run-total per-move-type counters: always reset, independent of the
    # cluster-move configuration above
    empty!(gc_params.move_stats)

    df = DataFrame(iter=Int[], omega=Float64[], energy=Float64[], num_particles=Int[])
    if observables !== nothing
        _validate_observables(observables, liveset)
        for (name, _) in observables
            df[!, name] = Float64[]
        end
    end
    culled = nothing

    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)

        if observables !== nothing || dead_point_callback !== nothing
            # Pre-sort with the step's own ordering (Ω = E − μN, descending,
            # cached keys through a stable sortperm, order-identical to the
            # by-comparator sort) and hold the walker the step will cull;
            # see `nested_sampling`.
            pre_keys = [_grand_potential(w, gc_params.chemical_potential)
                        for w in liveset.walkers]
            permute!(liveset.walkers, sortperm(pre_keys, rev=true))
            culled = liveset.walkers[1]
        end

        iter, omega, energy, n_par, liveset, gc_params = nested_sampling_step!(
            liveset, gc_params, mc_routine; ns_iteration=i)

        @debug "GC-NS step $i, iter: $iter, omega: $omega, energy: $energy, N: $n_par"

        if gc_params.fail_count >= gc_params.allowed_fail_count
            @warn "GC-NS: Failed $(gc_params.allowed_fail_count) times in a row."
            # Opt-in stall stop: return the partial ledger and the intact
            # live set instead of burning the remaining iteration budget
            # (the driver re-initializes its live set on entry, so a stalled
            # run cannot be chunked around from outside). The break leaves
            # fail_count at threshold so a caller can distinguish a stalled
            # return (the atomistic convention).
            stop_on_stall && break
            gc_params.fail_count = 0
        end

        if !(iter isa typeof(missing))
            if culled !== nothing
                omega == _grand_potential(culled, gc_params.chemical_potential) || error(
                    "GC-NS observable recording: dead-point/observable " *
                    "pairing lost (step culled a walker with Ω = $omega, but " *
                    "the pre-sorted worst walker had Ω = " *
                    "$(_grand_potential(culled, gc_params.chemical_potential)))")
            end
            if observables === nothing
                push!(df, (iter, omega.val, energy.val, n_par))
            else
                push!(df, (iter, omega.val, energy.val, n_par,
                           (Float64(f(culled.configuration)) for (_, f) in observables)...))
            end
            dead_point_callback === nothing || dead_point_callback(iter, culled)
        end

        if print_info && !(iter isa typeof(missing))
            @info "GC-NS iter: $(iter), Ω: $(omega), E: $(energy), N: $(n_par)"
        elseif print_info && iter isa typeof(missing)
            @info "GC-NS MC move failed, step: $(i)"
        end

        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end

    return df, liveset, gc_params
end

# ======================================================================
# Ideal-gas-referenced grand-canonical nested sampling
# ======================================================================

"""
    mutable struct IdealGasReferencedGCNSParameters <: SamplingParameters

Parameters for ideal-gas-referenced grand-canonical nested sampling on
lattice systems.

The prior is the non-interacting (ideal) lattice gas at reference fugacity
`z0`: a Bernoulli product measure in which each site is occupied independently
with probability `p0 = z0/(1 + z0)`, so a configuration with N particles
carries prior weight `z0^N` and the total prior mass on M sites is
`(1 + z0)^M`. Nested sampling culls on the energy E alone — the chemical
potential does not enter the sampler. Both μ and T are recovered in
post-processing (see `gc_thermodynamic_stats_ideal_ref` in `AnalysisTools`),
so a single run yields Ξ(μ, T) over a continuous temperature range and a
neighborhood of μ around `μ_ref(T) = k_B T ln z0`.

Setting `reference_fugacity = 1` makes the prior the uniform measure over all
2^M microstates — the same prior as the Ω-sorted construction in
`GrandCanonicalNestedSamplingParameters`, which instead bakes a single μ into
the sort quantity Ω = E − μN.

Unlike the Ω-sorted construction there is deliberately no `n_max` field: the
post-processing normalization `(1 + z0)^M` is the prior mass of the *full*
occupation range N ∈ {0, …, M}, and capping insertions would silently
truncate the prior support and bias Ξ.

# Fields
- `mc_steps::Int64`: MCMC steps per replacement walker.
- `reference_fugacity::Float64`: Reference fugacity z0 of the ideal-lattice-gas prior.
- `energy_perturbation::Float64`: Perturbation to break energy degeneracies.
  Required to be nonzero on lattices (degenerate levels stall the strict `<`
  ceiling and bias the evidence — enforced by the keyword constructor); must
  remain ≪ k_B·T at the lowest temperature targeted in post-processing, since
  perturbed energies are recorded. It must also stay *above* the float
  resolution of the largest energies on the ladder — keep
  `energy_perturbation / eps(E_max)` above ~`K²` so the `K` walkers draw
  distinct tie-breaking values on a degenerate plateau. The default `1e-12`
  satisfies that only for `E_max` up to `O(1) eV`; finite-J hard-core ladders
  reach `E_max = J·M·c/2`, so use `1e-9` there (safe through
  `E_max ≈ 10³ eV`, i.e. a few hundred sites at `J = 1 eV`) and scale it
  proportionally beyond — see the hard-core recipe in the
  `GenericLatticeHamiltonian` docstring.
- `random_seed::Int64`: Kept for parity with `GrandCanonicalNestedSamplingParameters`;
  **not currently consumed** by the NS loop — call `Random.seed!` before
  `ideal_gas_referenced_nested_sampling` for reproducible runs.
- `fail_count::Int64`: Consecutive failed replacements.
- `allowed_fail_count::Int64`: Maximum consecutive failures before warning.
- `cluster_p::Float64`: Current cluster growth probability (mutable runtime state).
- `cluster_accepted::Float64`: Accepted cluster moves in current adjustment window.
- `cluster_total::Float64`: Total cluster moves attempted in current adjustment window.
- `cluster_p_history::Vector{Float64}`: Trajectory of cluster_p after each adjustment.
- `cluster_accept_history::Vector{Float64}`: Acceptance rate at each adjustment.
- `cluster_adjust_iterations::Vector{Int}`: NS iteration index at each adjustment.
- `move_stats::Dict{Symbol,Int}`: Run-total per-move-type attempt/accept counters
  accumulated from every decorrelation walk (keys match the walk's `move_stats`
  NamedTuple; cleared once at run start, never window-reset).
"""
mutable struct IdealGasReferencedGCNSParameters <: SamplingParameters
    mc_steps::Int64
    reference_fugacity::Float64
    energy_perturbation::Float64
    random_seed::Int64
    fail_count::Int64
    allowed_fail_count::Int64
    cluster_p::Float64
    cluster_accepted::Float64
    cluster_total::Float64
    cluster_p_history::Vector{Float64}
    cluster_accept_history::Vector{Float64}
    cluster_adjust_iterations::Vector{Int}
    move_stats::Dict{Symbol,Int}
end

"""
    IdealGasReferencedGCNSParameters(;
        mc_steps=100, reference_fugacity=1.0, energy_perturbation=1e-12,
        random_seed=1234, fail_count=0, allowed_fail_count=10,
        cluster_p=0.3, cluster_accepted=0.0, cluster_total=0.0,
        cluster_p_history=Float64[], cluster_accept_history=Float64[],
        cluster_adjust_iterations=Int[], move_stats=Dict{Symbol,Int}())

Convenience constructor for `IdealGasReferencedGCNSParameters`.

`reference_fugacity` (z0) must be positive. Choose z0 near the target
fugacity range `exp(βμ)` of interest: post-run reweighting to a target (μ, T)
carries a factor `(exp(βμ)/z0)^N` per sample, and its effective sample size
degrades as `|βμ − ln z0|` grows (roughly beyond `1/√Var(N)`).

The `cluster_*` fields are mutable runtime state for adaptive cluster move
tuning, initialized from the static configuration on `MCGrandCanonicalMoves`
at the start of `ideal_gas_referenced_nested_sampling` when `clusters_freq > 0`.

`move_stats` holds run-total per-move-type attempt/accept counters accumulated
from every decorrelation walk (keys match the walk's `move_stats` NamedTuple;
never window-reset, cleared once at the start of each run).
"""
function IdealGasReferencedGCNSParameters(;
    mc_steps::Int64=100,
    reference_fugacity::Float64=1.0,
    energy_perturbation::Float64=1e-12,
    random_seed::Int64=1234,
    fail_count::Int64=0,
    allowed_fail_count::Int64=10,
    cluster_p::Float64=0.3,
    cluster_accepted::Float64=0.0,
    cluster_total::Float64=0.0,
    cluster_p_history::Vector{Float64}=Float64[],
    cluster_accept_history::Vector{Float64}=Float64[],
    cluster_adjust_iterations::Vector{Int}=Int[],
    move_stats::Dict{Symbol,Int}=Dict{Symbol,Int}(),
)
    if reference_fugacity <= 0.0
        throw(ArgumentError("reference_fugacity must be positive"))
    end
    if energy_perturbation == 0.0
        throw(ArgumentError("energy_perturbation must be nonzero: degenerate " *
            "lattice energy levels stall the strict < ceiling and bias the evidence"))
    end
    IdealGasReferencedGCNSParameters(
        mc_steps, reference_fugacity, energy_perturbation,
        random_seed, fail_count, allowed_fail_count,
        cluster_p, cluster_accepted, cluster_total,
        cluster_p_history, cluster_accept_history, cluster_adjust_iterations,
        move_stats,
    )
end

"""
    _init_ideal_gas_ref_walkers!(liveset::LatticeGasWalkers,
                                 params::IdealGasReferencedGCNSParameters)

Initialize walkers as exact i.i.d. draws from the ideal-lattice-gas prior:
each site occupied independently with probability `z0/(1 + z0)`.
"""
function _init_ideal_gas_ref_walkers!(liveset::LatticeGasWalkers,
                                      params::IdealGasReferencedGCNSParameters)
    h = liveset.hamiltonian
    z0 = params.reference_fugacity
    p0 = z0 / (1.0 + z0)
    for walker in liveset.walkers
        random_microstate!(walker.configuration; p=p0)
        assign_energy!(walker, h; perturb_energy=params.energy_perturbation)
        # Reset the iteration counter: df.iter feeds the ωᵢ prior-volume
        # weights, so a stale counter from a reused liveset corrupts them
        walker.iter = 0
    end
    return liveset
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers,
                          params::IdealGasReferencedGCNSParameters,
                          mc_routine::MCGrandCanonicalMoves;
                          ns_iteration::Int=0)

Perform one step of ideal-gas-referenced grand-canonical nested sampling.

Sorts walkers by energy E (not Ω — the chemical potential plays no role in
the sampler), removes the worst (highest E), clones a parent with E < E_worst,
and decorrelates the clone via grand-canonical MCMC that preserves the
`z0^N`-weighted prior below the energy ceiling.

# Returns
- `iter`: Iteration number (or `missing` if the step failed).
- `emax`: The E value of the removed walker (with units).
- `num_particles`: The N value of the removed walker.
- `liveset`: The updated liveset.
- `params`: The updated parameters.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers,
                               params::IdealGasReferencedGCNSParameters,
                               mc_routine::MCGrandCanonicalMoves;
                               ns_iteration::Int=0)
    ats = liveset.walkers
    h = liveset.hamiltonian
    n_walkers = length(ats)

    # Sort by energy (descending — worst first); the z0^N prior weighting is
    # maintained by the MC walk, not by the sort
    sort_by_energy!(liveset)

    iter::Union{Missing,Int} = missing
    worst = ats[1]
    emax_worst = worst.energy
    n_worst = sum(worst.configuration.components[1])

    emax_val = emax_worst.val  # unitless for the MC function
    eligible = [k for k in 2:n_walkers if ats[k].energy < emax_worst]
    if !isempty(eligible)
        parent_idx = rand(eligible)
    else
        parent_idx = rand(2:n_walkers)
    end
    to_walk = _clone_walker_shared_geometry(ats[parent_idx])

    # Decorrelate via GC MCMC: with mu = 0 the Ω ceiling reduces to an energy
    # ceiling, and z0 weights the insert/delete acceptance to preserve the
    # Bernoulli(z0/(1+z0)) prior
    accept, rate, to_walk, cl_accepted, cl_total, move_stats = MC_grand_canonical_walk!(
        params.mc_steps, to_walk, h, emax_val, 0.0;
        p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
        energy_perturb=params.energy_perturbation,
        z0=params.reference_fugacity,
        clusters_freq=mc_routine.clusters_freq,
        swaps_freq=mc_routine.swaps_freq,
        cluster_p=params.cluster_p,
        p_bias=mc_routine.p_bias,
        bias_predicate=mc_routine.bias_predicate,
        bias_shells=mc_routine.bias_shells,
        incremental=mc_routine.incremental)

    if accept
        push!(ats, to_walk)
        popfirst!(ats)
        update_iter!(liveset)
        params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        emax_worst = missing
        n_worst = missing
        params.fail_count += 1
    end

    # Accumulate cluster acceptance stats and adapt cluster_p
    _accumulate_cluster_stats!(params, mc_routine, cl_accepted, cl_total, ns_iteration)
    _accumulate_move_stats!(params, move_stats)

    return iter, emax_worst, n_worst, liveset, params
end

"""
    ideal_gas_referenced_nested_sampling(liveset::LatticeGasWalkers,
                                         params::IdealGasReferencedGCNSParameters,
                                         n_steps::Int64,
                                         mc_routine::MCGrandCanonicalMoves,
                                         save_strategy::DataSavingStrategy)

Run the ideal-gas-referenced grand-canonical nested sampling loop.

Walkers are initialized as exact i.i.d. draws from the ideal-lattice-gas
prior at reference fugacity z0 (each site occupied with probability
`z0/(1 + z0)`), then the loop iterates: remove the highest-energy walker,
record (E, N), replace with a clone decorrelated below the energy ceiling
under the `z0^N`-weighted prior. The chemical potential never enters the
run; a single run is post-processed to Ξ(μ, T) on an arbitrary (μ, T) grid
by `gc_thermodynamic_stats_ideal_ref` in `AnalysisTools`, which also needs
the surviving live walkers' (E, N) for the live-set tail closure —
extract them from the returned liveset.

# Arguments
- `liveset::LatticeGasWalkers`: The initial liveset (walkers will be re-initialized).
- `params::IdealGasReferencedGCNSParameters`: Parameters including the reference fugacity z0.
- `n_steps::Int64`: Number of NS iterations.
- `mc_routine::MCGrandCanonicalMoves`: The GC move routine (reused from the Ω-sorted construction).
- `save_strategy::DataSavingStrategy`: Strategy for periodic output.
- `observables`: Optional vector of `name::Symbol => callback` pairs
  recording per-dead-point observables (e.g. `order_parameter_c2x2`) as
  extra ledger columns; see [`nested_sampling`](@ref). The recorded columns
  feed the `observable_cols` machinery of
  `gc_thermodynamic_stats_ideal_ref`.
- `dead_point_callback`: Optional per-dead-point callback
  `(iter, walker) -> ...` invoked with the culled walker immediately after
  its ledger row is pushed; see [`nested_sampling`](@ref) for the contract.
  Configurations collected here support post-run re-evaluation under a
  second Hamiltonian, whose energies substitute for the `emax` column (and
  the live tail) in `gc_thermodynamic_stats_ideal_ref`: the prior-volume
  weights depend only on the iteration index, so the substituted column is
  a valid estimator of the second Hamiltonian's grand potential over the
  same ladder, with the returned Kish N_eff as its fidelity diagnostic.

- `stop_on_stall::Bool=false`: When true and `fail_count` reaches
  `allowed_fail_count`, warn once and return the partial ledger and the
  intact live set (`fail_count` stays at threshold); the default keeps the
  shipped warn-and-continue behavior byte-identically.

# Returns
- `df::DataFrame`: Columns `[:iter, :emax, :num_particles]`,
  plus one column per requested observable.
- `liveset::LatticeGasWalkers`: The final liveset (surviving walkers).
- `params::IdealGasReferencedGCNSParameters`: Updated parameters.
"""
function ideal_gas_referenced_nested_sampling(liveset::LatticeGasWalkers,
                                              params::IdealGasReferencedGCNSParameters,
                                              n_steps::Int64,
                                              mc_routine::MCGrandCanonicalMoves,
                                              save_strategy::DataSavingStrategy;
                                              observables::Union{Nothing,AbstractVector{<:Pair{Symbol,<:Any}}}=nothing,
                                              dead_point_callback::Union{Nothing,Function}=nothing,
                                              stop_on_stall::Bool=false)
    # Initialize walkers as i.i.d. draws from the Bernoulli(z0/(1+z0)) prior
    _init_ideal_gas_ref_walkers!(liveset, params)
    _warn_perturbation_scale(liveset, params.energy_perturbation)

    # Initialize cluster_p and reset counters from MCGrandCanonicalMoves if applicable
    if mc_routine.clusters_freq > 0
        params.cluster_p = mc_routine.initial_cluster_p
        params.cluster_accepted = 0.0
        params.cluster_total = 0.0
        empty!(params.cluster_p_history)
        empty!(params.cluster_accept_history)
        empty!(params.cluster_adjust_iterations)
    end
    # Run-total per-move-type counters: always reset, independent of the
    # cluster-move configuration above
    empty!(params.move_stats)

    df = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[])
    if observables !== nothing
        _validate_observables(observables, liveset)
        for (name, _) in observables
            df[!, name] = Float64[]
        end
    end
    culled = nothing

    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)

        if observables !== nothing || dead_point_callback !== nothing
            # Pre-sort with the step's own comparator (energy, descending)
            # and hold the walker the step will cull; see `nested_sampling`.
            sort_by_energy!(liveset)
            culled = liveset.walkers[1]
        end

        iter, emax, n_par, liveset, params = nested_sampling_step!(
            liveset, params, mc_routine; ns_iteration=i)

        @debug "IG-ref GC-NS step $i, iter: $iter, emax: $emax, N: $n_par"

        if params.fail_count >= params.allowed_fail_count
            @warn "IG-ref GC-NS: Failed $(params.allowed_fail_count) times in a row."
            # Opt-in stall stop (see grand_canonical_nested_sampling); the
            # break leaves fail_count at threshold
            stop_on_stall && break
            params.fail_count = 0
        end

        if !(iter isa typeof(missing))
            if culled !== nothing
                emax == culled.energy || error(
                    "IG-ref GC-NS observable recording: dead-point/observable " *
                    "pairing lost (step culled a walker with energy $emax, " *
                    "but the pre-sorted worst walker had $(culled.energy))")
            end
            if observables === nothing
                push!(df, (iter, emax.val, n_par))
            else
                push!(df, (iter, emax.val, n_par,
                           (Float64(f(culled.configuration)) for (_, f) in observables)...))
            end
            dead_point_callback === nothing || dead_point_callback(iter, culled)
        end

        if print_info && !(iter isa typeof(missing))
            @info "IG-ref GC-NS iter: $(iter), E: $(emax), N: $(n_par)"
        elseif print_info && iter isa typeof(missing)
            @info "IG-ref GC-NS MC move failed, step: $(i)"
        end

        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end

    return df, liveset, params
end
"""
    struct MCGalileanWalk <: MCRoutine

Move routine for canonical atomistic nested sampling through the Galilean
reflective walk kernel `MC_galilean_walk!`: all free particles move
collectively along straight lines, reflecting off the constraint boundary
through the energy gradient. Clone semantics (a random survivor is cloned and
decorrelated); three-dimensional only. `ns_params.mc_steps` sets the number of
trajectories per replacement and `ns_params.step_size` the 3N-space segment
length, adapted through the standard step-size machinery on the segment rate.

# Fields
- `n_refresh::Int`: Straight-line segments per trajectory (velocity redrawn
  between trajectories).
"""
struct MCGalileanWalk <: MCRoutine
    n_refresh::Int
    function MCGalileanWalk(; n_refresh::Int=8)
        n_refresh > 0 || throw(ArgumentError("n_refresh must be positive"))
        new(n_refresh)
    end
end

"""
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters,
                          mc_routine::MCGalileanWalk; ns_iteration::Int=0)

One canonical nested-sampling step decorrelated by the Galilean reflective walk:
the shape of the serial random-walk step (plateau-aware tie eviction included,
with refill walks running the reflective kernel below the plateau), with the
replacement clone drawn from the survivors and decorrelated through
`MC_galilean_walk!`. Returns the serial steps' five-value shape, so the driver
records `log_compression` as usual.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCGalileanWalk; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    pot = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    log_t::Union{Missing,Float64} = missing
    n_live = length(ats)
    n_tied = _tie_block_length(ats)
    if (n_tied >= 2 || ns_params.plateau_refill_target != 0) && n_tied < n_live
        if ns_params.plateau_refill_target == 0
            ns_params.plateau_refill_target = n_live
        end
        popfirst!(ats)
        update_iter!(liveset)
        iter = liveset.walkers[1].iter
        log_t = log((n_live - 1) / n_live)
        if n_tied == 1
            refill_fails = 0
            while length(ats) < ns_params.plateau_refill_target && refill_fails < ns_params.allowed_fail_count
                accept_r, _, at_r = MC_galilean_walk!(ns_params.mc_steps, deepcopy(rand(ats)), pot, emax;
                                                      step_size=ns_params.step_size,
                                                      n_refresh=mc_routine.n_refresh)
                if accept_r
                    push!(ats, at_r)
                else
                    refill_fails += 1
                end
            end
            length(ats) < ns_params.plateau_refill_target && @warn "Plateau refill left the live set at $(length(ats)) of $(ns_params.plateau_refill_target) walkers after $(refill_fails) failed decorrelation attempts; subsequent culls are charged with the actual live count."
            ns_params.plateau_refill_target = 0
        end
        return iter, emax, liveset, ns_params, log_t
    end
    to_walk = deepcopy(rand(ats[2:end]))
    accept, rate, at = MC_galilean_walk!(ns_params.mc_steps, to_walk, pot, emax;
                                         step_size=ns_params.step_size,
                                         n_refresh=mc_routine.n_refresh)
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
        log_t = log(n_live / (n_live + 1))
    else
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate; range=ns_params.accept_range)
    return iter, emax, liveset, ns_params, log_t
end

"""
    struct MCAtomGrandCanonicalMoves <: MCRoutine

Move routine for atomistic ideal-gas-referenced grand-canonical nested sampling:
single-atom displacements mixed with continuous-space insertions and deletions through
the atomistic `MC_grand_canonical_walk!` kernel, optionally cavity-biased through the
kernel's deterministic grid-cell channel (`p_bias`, `bias_radius`, `bias_grid`,
forwarded verbatim; the defaults keep the uniform channel and its random stream
bit-identically). Deliberately lean otherwise: the lattice `MCGrandCanonicalMoves`
carries cluster configuration that has no continuous-space counterpart, and silently
ignored fields on a mismatched routine are a known hazard class.

# Fields
- `p_move::Float64`: Probability of a displacement move.
- `p_insert::Float64`: Probability of an insertion move (deletion takes the remainder).
- `step_rate_source::Symbol`: Which acceptance rate drives the step-size
  adjustment. `:mixed` (the default) keeps the historical behavior, adapting on
  the kernel's combined rate over all three channels; `:move` adapts on the
  displacement-only rate from the walk's per-channel counters, so a collapsing
  insertion or deletion acceptance in a dense interacting system cannot drag the
  displacement step size down against a healthy displacement channel (adjustment
  is skipped when a walk attempted no displacement).
- `mc_steps_per_particle::Float64`: Extra kernel steps per parent particle: a
  decorrelation walk for a parent at particle count N runs
  `mc_steps + round(Int, mc_steps_per_particle * N)` steps, restoring per-sweep
  decorrelation as the count grows. The default 0.0 is draw-count identical to
  the historical fixed-length walk.
"""
struct MCAtomGrandCanonicalMoves <: MCRoutine
    p_move::Float64
    p_insert::Float64
    step_rate_source::Symbol
    mc_steps_per_particle::Float64
    galilean_steps::Int
    galilean_n_refresh::Int
    galilean_step_size::Float64
    p_bias::Float64
    bias_radius::Float64
    bias_grid::Int
    function MCAtomGrandCanonicalMoves(; p_move::Float64=0.5, p_insert::Float64=0.25,
                                       step_rate_source::Symbol=:mixed,
                                       mc_steps_per_particle::Float64=0.0,
                                       galilean_steps::Int=0,
                                       galilean_n_refresh::Int=8,
                                       galilean_step_size::Float64=0.5,
                                       p_bias::Float64=0.0,
                                       bias_radius::Float64=0.0,
                                       bias_grid::Int=0)
        if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
            throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
        end
        if !(step_rate_source in (:mixed, :move))
            throw(ArgumentError("step_rate_source must be :mixed or :move"))
        end
        if !(mc_steps_per_particle >= 0.0)
            throw(ArgumentError("mc_steps_per_particle must be non-negative"))
        end
        if galilean_steps < 0 || galilean_n_refresh <= 0 || galilean_step_size <= 0.0
            throw(ArgumentError("galilean_steps must be non-negative and galilean_n_refresh, galilean_step_size positive"))
        end
        if p_bias < 0.0 || p_bias > 1.0
            throw(ArgumentError("p_bias must lie in [0, 1]"))
        end
        if p_bias > 0.0 && (bias_radius <= 0.0 || bias_grid < 1)
            throw(ArgumentError("a biased insertion channel (p_bias > 0) requires bias_radius > 0 and bias_grid >= 1"))
        end
        if p_bias == 1.0
            @warn "p_bias = 1: any deletion whose vacated cell is not a post-deletion cavity cell auto-rejects, which can freeze dense states; a mixed channel (p_bias < 1) keeps the chain irreducible." maxlog=1
        end
        new(p_move, p_insert, step_rate_source, mc_steps_per_particle,
            galilean_steps, galilean_n_refresh, galilean_step_size,
            p_bias, bias_radius, bias_grid)
    end
end

# Number of kernel steps for a decorrelation walk whose parent carries n_par
# particles (see `MCAtomGrandCanonicalMoves.mc_steps_per_particle`)
function _gc_walk_length(params::SamplingParameters, mc_routine::MCAtomGrandCanonicalMoves, n_par::Int)
    return params.mc_steps + round(Int, mc_routine.mc_steps_per_particle * n_par)
end

"""
    mutable struct AtomisticIGRefGCNSParameters <: SamplingParameters

Parameters for atomistic energy-sorted ideal-gas-referenced grand-canonical nested
sampling. The sampler is athermal: the chemical potential and the temperature never
enter a run; both enter only in the post-run reduction to Ξ(μ, T).

# Fields
- `mc_steps::Int64`: MCMC steps per replacement walker.
- `reference_activity::typeof(1.0u"Å^-3")`: The reference activity z0 (dimension inverse
  volume). The driver folds the dimensionless z0V once from the walkers' shared
  orthorhombic cell.
- `species::Symbol`: Chemical identity of inserted particles, passed explicitly because
  a walker's configuration may be empty.
- `initial_step_size::Float64`: Initial displacement step size (Angstrom).
- `step_size::Float64`: Current displacement step size (mutable runtime state).
- `step_size_lo::Float64`: Lower bound for the step-size adjustment.
- `step_size_up::Float64`: Upper bound for the step-size adjustment.
- `accept_range::Tuple{Float64,Float64}`: Acceptance-rate window forwarded to the
  step-size adjustment.
- `fail_count::Int64`: Consecutive failed replacements.
- `allowed_fail_count::Int64`: Maximum consecutive failures before the driver's stall
  contract fires (see `ideal_gas_referenced_nested_sampling`). The default (100,
  matching the canonical atomistic loop) suits interacting descents, where a run of
  failed replacements is an ordinary sampling hiccup rather than a terminal state.
- `plateau_refill_target::Int64`: Live-set size to restore after a plateau eviction
  block (mutable runtime state; reset at driver entry).
- `refill_fail_budget::Int64`: Failure budget for one plateau-refill loop. The default
  0 keeps the historical behavior of charging refill failures against
  `allowed_fail_count` cumulatively; a positive value gives each refill loop its own
  budget, so recovering from one difficult plateau block cannot permanently shrink
  the live set within the stall budget.
- `move_stats::Dict{Symbol,Int}`: Run-total per-move-type attempt/accept counters
  accumulated from every ordinary-cull decorrelation walk (cleared once at run
  start; plateau-refill walks are not accumulated).

Three fields carried by sibling parameter structs are deliberately absent. No `n_max`:
the reference measure has unbounded particle-number support, and a cap silently
truncates its mass and biases the evidence (the kernel's keyword exists for bounded
constructions only; this sampler never sets it). No `energy_perturbation`: exact energy
ties are handled by the plateau-aware eviction machinery, and perturbing recorded
energies would break the exactly-zero closed-form checks this construction is validated
against. No `random_seed`: the lattice ideal-gas-referenced sibling documents its field
as not consumed by the sampling loop; seed the global RNG before the run instead.
"""
mutable struct AtomisticIGRefGCNSParameters <: SamplingParameters
    mc_steps::Int64
    reference_activity::typeof(1.0u"Å^-3")
    species::Symbol
    initial_step_size::Float64
    step_size::Float64
    step_size_lo::Float64
    step_size_up::Float64
    accept_range::Tuple{Float64,Float64}
    fail_count::Int64
    allowed_fail_count::Int64
    plateau_refill_target::Int64
    refill_fail_budget::Int64
    move_stats::Dict{Symbol,Int}
end

"""
    AtomisticIGRefGCNSParameters(;
        mc_steps=100, reference_activity=0.01u"Å^-3", species=:H,
        initial_step_size=0.5, step_size=0.5, step_size_lo=0.01, step_size_up=2.0,
        accept_range=(0.25, 0.75), fail_count=0, allowed_fail_count=100,
        plateau_refill_target=0, refill_fail_budget=0, move_stats=Dict{Symbol,Int}())

Convenience constructor for `AtomisticIGRefGCNSParameters`. `reference_activity` (z0)
must be positive; choose it so z0V sits near the particle-number range of interest, as
post-run reweighting to a target (μ, T) carries a factor `(z/z0)^N` per sample whose
effective sample size degrades away from the reference.
"""
function AtomisticIGRefGCNSParameters(;
    mc_steps::Int64=100,
    reference_activity::typeof(1.0u"Å^-3")=0.01u"Å^-3",
    species::Symbol=:H,
    initial_step_size::Float64=0.5,
    step_size::Float64=0.5,
    step_size_lo::Float64=0.01,
    step_size_up::Float64=2.0,
    accept_range::Tuple{Float64,Float64}=(0.25, 0.75),
    fail_count::Int64=0,
    allowed_fail_count::Int64=100,
    plateau_refill_target::Int64=0,
    refill_fail_budget::Int64=0,
    move_stats::Dict{Symbol,Int}=Dict{Symbol,Int}(),
)
    if reference_activity <= 0.0u"Å^-3"
        throw(ArgumentError("reference_activity must be positive"))
    end
    if refill_fail_budget < 0
        throw(ArgumentError("refill_fail_budget must be non-negative (0 charges refill failures against allowed_fail_count)"))
    end
    AtomisticIGRefGCNSParameters(
        mc_steps, reference_activity, species,
        initial_step_size, step_size, step_size_lo, step_size_up, accept_range,
        fail_count, allowed_fail_count, plateau_refill_target, refill_fail_budget, move_stats,
    )
end

# Per-iteration acceptance-ledger columns (order matters: it is the ledger
# schema order when `record_move_rates=true`)
const _MOVE_RATE_COLUMNS = (:move_attempted, :move_accepted, :insert_attempted,
                            :insert_accepted, :delete_attempted, :delete_accepted)

_move_stats_snapshot(params::SamplingParameters) =
    Tuple(get(params.move_stats, k, 0) for k in _MOVE_RATE_COLUMNS)

"""
    _warn_min_image_cutoff(pot, config)

Warn once when a potential's finite interaction range exceeds half the smallest
cell edge of `config`: minimum-image truncation is anisotropic in that regime.
An infinite range never warns (an untruncated minimum-image Hamiltonian is a
deliberate model choice), and potentials of undeterminable range are skipped.
"""
function _warn_min_image_cutoff(pot, config)
    rng = AbstractPotentials._max_interaction_range(pot)
    rng === missing && return nothing
    isfinite(ustrip(rng)) || return nothing
    cellv = cell_vectors(config)
    L_min = min(ustrip(u"Å", cellv[1][1]), ustrip(u"Å", cellv[2][2]), ustrip(u"Å", cellv[3][3]))
    if ustrip(u"Å", rng) > L_min / 2
        @warn "The potential's interaction range ($(round(ustrip(u"Å", rng); sigdigits=4)) Å) exceeds half the smallest cell edge ($(round(L_min / 2; sigdigits=4)) Å): minimum-image truncation is anisotropic in this regime."
    end
    return nothing
end

"""
    _atomistic_igref_z0V(liveset::AtomWalkers, params::AtomisticIGRefGCNSParameters)

Validate the liveset for the atomistic ideal-gas-referenced construction (single
unfrozen component per walker, one shared orthorhombic cell) and return the
dimensionless product of the reference activity and the cell volume.
"""
function _atomistic_igref_z0V(liveset::AtomWalkers, params::AtomisticIGRefGCNSParameters)
    isempty(liveset.walkers) && throw(ArgumentError("the liveset carries no walkers"))
    cellv = cell_vectors(liveset.walkers[1].configuration)
    for i in 1:3, j in 1:3
        if i != j && !iszero(ustrip(cellv[i][j]))
            throw(ArgumentError("the atomistic ideal-gas-referenced construction assumes an orthorhombic cell (consistent with pbc_dist); found a nonzero off-diagonal cell component"))
        end
    end
    for walker in liveset.walkers
        walker isa AtomWalker{1} || throw(ArgumentError("every walker must be a single-component AtomWalker{1}"))
        any(walker.frozen) && throw(ArgumentError("the atomistic ideal-gas-referenced construction requires unfrozen walkers"))
        cell_vectors(walker.configuration) == cellv || throw(ArgumentError("all walkers must share one cell"))
    end
    V = cellv[1][1] * cellv[2][2] * cellv[3][3]
    return ustrip(Unitful.NoUnits, params.reference_activity * V)
end

"""
    _init_atomistic_igref_walkers!(liveset::AtomWalkers,
                                   params::AtomisticIGRefGCNSParameters,
                                   z0V::Float64)

Initialize walkers as exact i.i.d. draws from the continuous ideal-gas prior at the
reference activity: each walker's particle count is Poisson(z0V) and its positions are
uniform in the cell.
"""
function _init_atomistic_igref_walkers!(liveset::AtomWalkers,
                                        params::AtomisticIGRefGCNSParameters,
                                        z0V::Float64)
    pot = liveset.potential
    cellv = cell_vectors(liveset.walkers[1].configuration)
    box = (cellv[1][1], cellv[2][2], cellv[3][3])
    for walker in liveset.walkers
        while walker.list_num_par[1] > 0
            remove_particle!(walker, walker.list_num_par[1])
        end
        n = rand(Poisson(z0V))
        for _ in 1:n
            pos = SVector(rand() * box[1], rand() * box[2], rand() * box[3])
            insert_particle!(walker, pos, params.species)
        end
        AbstractLiveSets.assign_frozen_energy!(walker, pot)
        assign_energy!(walker, pot)
        # Reset the iteration counter: a stale counter from a reused liveset
        # corrupts the prior-volume bookkeeping
        walker.iter = 0
    end
    return liveset
end

"""
    nested_sampling_step!(liveset::AtomWalkers,
                          params::AtomisticIGRefGCNSParameters,
                          mc_routine::MCAtomGrandCanonicalMoves;
                          ns_iteration::Int=0, z0V::Union{Nothing,Float64}=nothing)

Perform one step of atomistic energy-sorted ideal-gas-referenced grand-canonical nested
sampling: sort by energy (descending), handle exact ceiling ties with the plateau-aware
eviction machinery of the canonical atomistic steps (refill walks run the grand-canonical
kernel, so refilled clones may re-enter at a different particle count; the compression
bookkeeping counts walkers, not particle-number sectors), otherwise cull the worst walker
and decorrelate a strictly-below-ceiling clone through `MC_grand_canonical_walk!` under
the z0-weighted prior. Accepted clones (cull and refill alike) have their energies
re-anchored by a from-scratch `interacting_energy` recompute and are re-checked against
the ceiling before entering the live set: the kernel's incremental bookkeeping leaves
rounding dust that would otherwise fragment bit-exact energy plateaus into artificial
sub-levels and admit clones whose true energy sits on, not below, the ceiling.

# Returns
- `iter`: Iteration number (or `missing` if the step failed).
- `emax`: The energy of the culled walker (with units; `missing` on failure).
- `num_particles`: The particle count of the culled walker (`missing` on failure).
- `liveset`: The updated liveset.
- `params`: The updated parameters.
- `log_t`: The log prior-volume compression charged for this cull (`missing` on failure).
"""
function nested_sampling_step!(liveset::AtomWalkers,
                               params::AtomisticIGRefGCNSParameters,
                               mc_routine::MCAtomGrandCanonicalMoves;
                               ns_iteration::Int=0,
                               z0V::Union{Nothing,Float64}=nothing)
    z0V === nothing && (z0V = _atomistic_igref_z0V(liveset, params))
    sort_by_energy!(liveset)
    ats = liveset.walkers
    pot = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = ats[1].energy
    num_particles::Union{Missing,Int} = ats[1].list_num_par[1]
    log_t::Union{Missing,Float64} = missing
    n_live = length(ats)
    n_tied = _tie_block_length(ats)
    if (n_tied >= 2 || params.plateau_refill_target != 0) && n_tied < n_live
        # Plateau block: evict the worst tied walker without replacement
        # (Fowlie-Handley-Su), charging (n_live - 1)/n_live with the shrinking
        # live count; see the canonical atomistic step for the full contract.
        if params.plateau_refill_target == 0
            params.plateau_refill_target = n_live
        end
        popfirst!(ats)
        update_iter!(liveset)
        iter = liveset.walkers[1].iter
        log_t = log((n_live - 1) / n_live)
        if n_tied == 1
            # The plateau is exhausted: refill by cloning survivors and
            # decorrelating strictly below the plateau energy with the
            # grand-canonical kernel (volume-neutral; no ledger row).
            refill_fails = 0
            refill_budget = params.refill_fail_budget > 0 ? params.refill_fail_budget : params.allowed_fail_count
            while length(ats) < params.plateau_refill_target && refill_fails < refill_budget
                at_r = deepcopy(rand(ats))
                accept_r, _, at_r, _ = MC_grand_canonical_walk!(
                    _gc_walk_length(params, mc_routine, at_r.list_num_par[1]), at_r, pot, emax;
                    z0V=z0V, species=params.species,
                    p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
                    step_size=params.step_size,
                    p_bias=mc_routine.p_bias, bias_radius=mc_routine.bias_radius,
                    bias_grid=mc_routine.bias_grid)
                if accept_r && mc_routine.galilean_steps > 0 && at_r.list_num_par[1] > 0
                    # Optional reflective decorrelation burst at the clone's
                    # current particle count; measure-preserving at fixed N, so
                    # the composition preserves the constrained reference
                    MC_galilean_walk!(mc_routine.galilean_steps, at_r, pot, emax;
                                      step_size=mc_routine.galilean_step_size,
                                      n_refresh=mc_routine.galilean_n_refresh)
                end
                if accept_r
                    # Re-anchor the clone's energy from scratch: the kernel's
                    # incremental bookkeeping leaves dust that would fragment
                    # bit-exact plateaus and admit sub-ceiling dust crossings
                    at_r.energy = interacting_energy(at_r.configuration, pot,
                        at_r.list_num_par, at_r.frozen) + at_r.energy_frozen_part
                    accept_r = at_r.energy < emax
                end
                if accept_r
                    push!(ats, at_r)
                else
                    refill_fails += 1
                end
            end
            length(ats) < params.plateau_refill_target && @warn "Plateau refill left the live set at $(length(ats)) of $(params.plateau_refill_target) walkers after $(refill_fails) failed decorrelation attempts; subsequent culls are charged with the actual live count."
            params.plateau_refill_target = 0
        end
        # Step-size adaptation and move-stats accumulation stay off inside
        # plateau blocks by design: refill walks run under a fixed plateau
        # ceiling, and their acceptance is not representative of ordinary culls
        return iter, emax, num_particles, liveset, params, log_t
    end
    # Ordinary cull: strictly-below-ceiling parent selection, mirroring the
    # lattice ideal-gas-referenced step
    eligible = [k for k in 2:n_live if ats[k].energy < emax]
    parent_idx = isempty(eligible) ? rand(2:n_live) : rand(eligible)
    to_walk = deepcopy(ats[parent_idx])
    accept, rate, to_walk, move_stats = MC_grand_canonical_walk!(
        _gc_walk_length(params, mc_routine, to_walk.list_num_par[1]), to_walk, pot, emax;
        z0V=z0V, species=params.species,
        p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
        step_size=params.step_size,
        p_bias=mc_routine.p_bias, bias_radius=mc_routine.bias_radius,
        bias_grid=mc_routine.bias_grid)
    if accept && mc_routine.galilean_steps > 0 && to_walk.list_num_par[1] > 0
        # Optional reflective decorrelation burst (see the refill branch)
        MC_galilean_walk!(mc_routine.galilean_steps, to_walk, pot, emax;
                          step_size=mc_routine.galilean_step_size,
                          n_refresh=mc_routine.galilean_n_refresh)
    end
    if accept
        # Re-anchor the clone's energy from scratch (see the refill branch)
        to_walk.energy = interacting_energy(to_walk.configuration, pot,
            to_walk.list_num_par, to_walk.frozen) + to_walk.energy_frozen_part
        accept = to_walk.energy < emax
    end
    if accept
        push!(ats, to_walk)
        popfirst!(ats)
        update_iter!(liveset)
        params.fail_count = 0
        iter = liveset.walkers[1].iter
        log_t = log(n_live / (n_live + 1))
    else
        emax = missing
        num_particles = missing
        params.fail_count += 1
    end
    _accumulate_move_stats!(params, move_stats)
    if mc_routine.step_rate_source === :move
        # Displacement-only adaptation: skip when the walk attempted no
        # displacement (nothing to adapt on)
        if move_stats.move_attempted > 0
            adjust_step_size(params, move_stats.move_accepted / move_stats.move_attempted;
                             range=params.accept_range)
        end
    else
        adjust_step_size(params, rate; range=params.accept_range)
    end
    return iter, emax, num_particles, liveset, params, log_t
end

"""
    ideal_gas_referenced_nested_sampling(liveset::AtomWalkers,
                                         params::AtomisticIGRefGCNSParameters,
                                         n_steps::Int64,
                                         mc_routine::MCAtomGrandCanonicalMoves,
                                         save_strategy::DataSavingStrategy;
                                         observables=nothing,
                                         dead_point_callback=nothing,
                                         stop_on_stall::Bool=true)

Run the atomistic energy-sorted ideal-gas-referenced grand-canonical nested sampling
loop. Walkers are initialized as exact i.i.d. draws from the continuous ideal-gas prior
at the reference activity (particle counts Poisson(z0V), positions uniform in the cell),
then the loop iterates: remove the highest-energy walker, record (E, N, log-compression),
replace with a strictly-below-ceiling clone decorrelated by the grand-canonical kernel
under the z0-weighted prior. Exact ceiling ties are handled by the plateau-aware
eviction machinery. Neither the chemical potential nor the temperature enters the run.

Stall contract: when `params.fail_count` reaches `params.allowed_fail_count`, the driver
warns once and, with `stop_on_stall=true` (the default), returns the partial ledger and
the intact live set instead of consuming the remaining iteration budget on proposals
that cannot be accepted. A fully degenerate live set (every configuration at exactly one
energy, the zero-interaction limit) is a legitimate terminal state whose entire prior
mass sits in the live set; the subsequent ledger reduction is then exact. With
`stop_on_stall=false` the driver keeps the warn-and-continue contract of the sibling
loops (the failure counter is reset, the step size is reset to
`params.initial_step_size` following the canonical loop's contract, and the loop
continues).

With `record_move_rates=true` the ledger gains the per-iteration acceptance columns
`[:move_attempted, :move_accepted, :insert_attempted, :insert_accepted,
:delete_attempted, :delete_accepted]`, recorded as the change in the run-total
counters since the previous recorded row: the production diagnostic for how the
per-channel acceptances evolve with descent depth. A recorded row therefore also
carries the attempts of any failed replacements since the previous row;
plateau-eviction rows run no decorrelation walk of their own, and refill-walk
counters are deliberately not accumulated. The recorded columns sum exactly to
the `move_stats` run totals accumulated through the last recorded row. The names are reserved ledger columns, rejected as
observable names.

Observable callbacks must handle empty configurations: under this construction a
walker's particle count may be zero.

# Returns
- `df::DataFrame`: Columns `[:iter, :emax, :num_particles, :log_compression]`, plus the
  acceptance columns when requested, plus one column per requested observable.
  Zero-accept runs return the schema with no rows.
- `liveset::AtomWalkers`: The final liveset (surviving walkers).
- `params::AtomisticIGRefGCNSParameters`: Updated parameters.
"""
function ideal_gas_referenced_nested_sampling(liveset::AtomWalkers,
                                              params::AtomisticIGRefGCNSParameters,
                                              n_steps::Int64,
                                              mc_routine::MCAtomGrandCanonicalMoves,
                                              save_strategy::DataSavingStrategy;
                                              observables::Union{Nothing,AbstractVector{<:Pair{Symbol,<:Any}}}=nothing,
                                              dead_point_callback::Union{Nothing,Function}=nothing,
                                              stop_on_stall::Bool=true,
                                              record_move_rates::Bool=false)
    z0V = _atomistic_igref_z0V(liveset, params)
    _warn_min_image_cutoff(liveset.potential, liveset.walkers[1].configuration)
    _init_atomistic_igref_walkers!(liveset, params, z0V)

    # Parameter-reuse hazards: runtime state never leaks between runs
    params.plateau_refill_target = 0
    params.fail_count = 0
    empty!(params.move_stats)

    df = DataFrame(iter=Int[], emax=Float64[], num_particles=Int[], log_compression=Float64[])
    if record_move_rates
        for name in _MOVE_RATE_COLUMNS
            df[!, name] = Int[]
        end
    end
    if observables !== nothing
        _validate_observables(observables, liveset)
        for (name, _) in observables
            df[!, name] = Float64[]
        end
    end
    culled = nothing
    empty_frame_warned = false
    rate_prev = record_move_rates ? _move_stats_snapshot(params) : nothing

    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        # The save layer cannot represent a zero-atom frame, and this
        # construction produces empty walkers by design (the reference measure
        # puts mass e^{-z0V} on N = 0): skip those writes with one warning
        if length(liveset.walkers[1].configuration) > 0
            write_walker_every_n(liveset.walkers[1], i, save_strategy)
        elseif !empty_frame_warned
            @warn "Skipping trajectory/live-set writes that would contain a zero-atom frame: the save layer cannot represent one."
            empty_frame_warned = true
        end

        if observables !== nothing || dead_point_callback !== nothing
            # Pre-sort with the step's own comparator (energy, descending)
            # and hold the walker the step will cull; see `nested_sampling`.
            sort_by_energy!(liveset)
            culled = liveset.walkers[1]
        end

        iter, emax, n_par, liveset, params, log_t = nested_sampling_step!(
            liveset, params, mc_routine; ns_iteration=i, z0V=z0V)

        @debug "Atomistic IG-ref GC-NS step $i, iter: $iter, emax: $emax, N: $n_par"

        if params.fail_count >= params.allowed_fail_count
            if stop_on_stall
                @warn "Atomistic IG-ref GC-NS: $(params.allowed_fail_count) consecutive failed replacements; returning the partial ledger and the intact live set (stop_on_stall=true)."
                break
            else
                @warn "Atomistic IG-ref GC-NS: Failed $(params.allowed_fail_count) times in a row. Reset step size!"
                params.fail_count = 0
                params.step_size = params.initial_step_size
            end
        end

        if !(iter isa typeof(missing))
            if culled !== nothing
                emax == culled.energy || error(
                    "Atomistic IG-ref GC-NS observable recording: dead-point/observable " *
                    "pairing lost (step culled a walker with energy $emax, " *
                    "but the pre-sorted worst walker had $(culled.energy))")
            end
            if record_move_rates
                snap = _move_stats_snapshot(params)
                rate_row = snap .- rate_prev
                rate_prev = snap
            else
                rate_row = ()
            end
            if observables === nothing
                push!(df, (iter, emax.val, n_par, log_t, rate_row...))
            else
                push!(df, (iter, emax.val, n_par, log_t, rate_row...,
                           (Float64(f(culled.configuration)) for (_, f) in observables)...))
            end
            dead_point_callback === nothing || dead_point_callback(iter, culled)
        end

        if print_info && !(iter isa typeof(missing))
            @info "Atomistic IG-ref GC-NS iter: $(iter), E: $(emax), N: $(n_par)"
        elseif print_info && iter isa typeof(missing)
            @info "Atomistic IG-ref GC-NS MC move failed, step: $(i)"
        end

        write_df_every_n(df, i, save_strategy)
        if all(w -> length(w.configuration) > 0, liveset.walkers)
            write_ls_every_n(liveset, i, save_strategy)
        elseif !empty_frame_warned
            @warn "Skipping trajectory/live-set writes that would contain a zero-atom frame: the save layer cannot represent one."
            empty_frame_warned = true
        end
    end

    return df, liveset, params
end
