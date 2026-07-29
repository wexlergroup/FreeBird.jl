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
            )
    NestedSamplingParameters(mc_steps, initial_step_size, step_size, step_size_lo, step_size_up, accept_range, fail_count, allowed_fail_count, energy_perturbation, random_seed, cluster_p, cluster_accepted, cluster_total, cluster_p_history, cluster_accept_history, cluster_adjust_iterations)
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

When `clusters_freq == 0` (the default), the fixed-N branch uses only local swaps
(`lattice_random_walk!`), preserving backward compatibility with existing scripts.
When `clusters_freq > 0`, the fixed-N branch mixes geometric cluster moves with
local swaps according to the `clusters_freq:swaps_freq` ratio.
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
    function MCGrandCanonicalMoves(;
            p_move::Float64=0.5,
            p_insert::Float64=0.25,
            clusters_freq::Int=0,
            swaps_freq::Int=1,
            initial_cluster_p::Float64=0.3,
            target_cluster_accept::Float64=0.3,
            cluster_adjust_interval::Int=50,
            cluster_p_floor::Float64=0.01,
            cluster_p_ceiling::Float64=1.0)
        if p_move < 0.0 || p_insert < 0.0 || p_move + p_insert > 1.0
            throw(ArgumentError("p_move and p_insert must satisfy 0 <= p_move + p_insert <= 1"))
        end
        new(p_move, p_insert, clusters_freq, swaps_freq,
            initial_cluster_p, target_cluster_accept, cluster_adjust_interval,
            cluster_p_floor, cluster_p_ceiling)
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
end

"""
    GrandCanonicalNestedSamplingParameters(;
        mc_steps=100, chemical_potential=0.0, energy_perturbation=1e-12,
        random_seed=1234, fail_count=0, allowed_fail_count=10,
        init_occupation_p=0.5, n_max=typemax(Int64),
        cluster_p=0.3, cluster_accepted=0.0, cluster_total=0.0,
        cluster_p_history=Float64[], cluster_accept_history=Float64[],
        cluster_adjust_iterations=Int[])

Convenience constructor for `GrandCanonicalNestedSamplingParameters`.

The `n_max` parameter sets an upper bound on the number of particles per walker.
Insertions are rejected when N ≥ n_max. Default is `typemax(Int64)` (no cap).

The `cluster_*` fields are mutable runtime state for adaptive cluster move tuning.
They are initialized from the static configuration on `MCGrandCanonicalMoves` at
the start of `grand_canonical_nested_sampling` when `clusters_freq > 0`.
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
)
    GrandCanonicalNestedSamplingParameters(
        mc_steps, chemical_potential, energy_perturbation,
        random_seed, fail_count, allowed_fail_count,
        init_occupation_p, n_max,
        cluster_p, cluster_accepted, cluster_total,
        cluster_p_history, cluster_accept_history, cluster_adjust_iterations,
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
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine)

Perform a single step of the nested sampling algorithm using the Monte Carlo random walk routine.

# Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRoutine`: The Monte Carlo routine for generating new samples. See [`MCRoutine`](@ref).

# Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
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
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate)
    return iter, emax, liveset, ns_params
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
    elseif length(mc_routine.dims) == 2
        random_walk_function = MC_random_walk_2D!
    else
        error("Unsupported dimensions: $(mc_routine.dims)")
    end

    mc_steps_per_worker = ceil(Int, ns_params.mc_steps / nworkers()) # distribute the total MC steps among workers

    walking = [remotecall(random_walk_function, workers()[i], mc_steps_per_worker, to_walk, lj, ns_params.step_size, emax[mc_routine.n_cull]) for (i,to_walk) in enumerate(to_walks)]
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

    adjust_step_size(ns_params, rate)
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
    elseif length(mc_routine.dims) == 2
        random_walk_function = MC_random_walk_2D!
    else
        error("Unsupported dimensions: $(mc_routine.dims)")
    end


    walking = [remotecall(random_walk_function, workers()[i], ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax[end]) for (i,to_walk) in enumerate(to_walks)]
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

    adjust_step_size(ns_params, rate)
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
    elseif length(mc_routine.dims) == 2
        random_walk_function = MC_random_walk_2D!
    else
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

    adjust_step_size(ns_params, rate)
    return iter, emax[end], liveset, ns_params
end

function nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine; ns_iteration::Int=0)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
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
    else
        # @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate)
    return iter, emax, liveset, ns_params
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
    adjust_step_size(ns_params, rate)

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

    adjust_step_size(ns_params, rate)
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

# Returns
- `df`: A DataFrame containing the iteration number and maximum energy for each step.
- `liveset`: The updated set of walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling(liveset::AbstractLiveSet, 
                                ns_params::NestedSamplingParameters, 
                                n_steps::Int64, 
                                mc_routine::MCRoutine,
                                save_strategy::DataSavingStrategy)
    # Initialize cluster_p and reset counters from MCMixedMoves if applicable
    if mc_routine isa MCMixedMoves && mc_routine.clusters_freq > 0
        ns_params.cluster_p = mc_routine.initial_cluster_p
        ns_params.cluster_accepted = 0.0
        ns_params.cluster_total = 0.0
        empty!(ns_params.cluster_p_history)
        empty!(ns_params.cluster_accept_history)
        empty!(ns_params.cluster_adjust_iterations)
    end
    df = DataFrame(iter=Int[], emax=Float64[])
    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)
        iter, emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine; ns_iteration=i)
        @debug "n_step $i, iter: $iter, emax: $emax"
        if ns_params.fail_count >= ns_params.allowed_fail_count
            @warn "Failed to accept MC move $(ns_params.allowed_fail_count) times in a row. Reset step size!"
            ns_params.fail_count = 0
            ns_params.step_size = ns_params.initial_step_size
        end
        if !(iter isa typeof(missing))
            push!(df, (iter, emax.val))
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

    # Sort by grand potential (descending — worst first)
    sort!(ats, by = w -> _grand_potential(w, mu), rev=true)

    iter::Union{Missing,Int} = missing
    worst = ats[1]
    omega_worst = _grand_potential(worst, mu)
    energy_worst = worst.energy
    n_worst = sum(worst.configuration.components[1])

    # Select parent: prefer walkers strictly below omega_worst
    omega_max_val = omega_worst.val  # unitless for the MC function
    eligible = [k for k in 2:n_walkers if _grand_potential(ats[k], mu) < omega_worst]
    if !isempty(eligible)
        parent_idx = rand(eligible)
    else
        parent_idx = rand(2:n_walkers)
    end
    to_walk = deepcopy(ats[parent_idx])

    # Decorrelate via GC MCMC
    accept, rate, to_walk, cl_accepted, cl_total = MC_grand_canonical_walk!(
        gc_params.mc_steps, to_walk, h, omega_max_val, mu;
        p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
        energy_perturb=gc_params.energy_perturbation,
        n_max=gc_params.n_max,
        clusters_freq=mc_routine.clusters_freq,
        swaps_freq=mc_routine.swaps_freq,
        cluster_p=gc_params.cluster_p)

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

# Returns
- `df::DataFrame`: Columns `[:iter, :omega, :energy, :num_particles]`.
- `liveset::LatticeGasWalkers`: The final liveset (surviving walkers).
- `gc_params::GrandCanonicalNestedSamplingParameters`: Updated parameters.
"""
function grand_canonical_nested_sampling(liveset::LatticeGasWalkers,
                                         gc_params::GrandCanonicalNestedSamplingParameters,
                                         n_steps::Int64,
                                         mc_routine::MCGrandCanonicalMoves,
                                         save_strategy::DataSavingStrategy)
    # Initialize walkers with random microstates
    _init_gc_walkers!(liveset, gc_params)

    # Initialize cluster_p and reset counters from MCGrandCanonicalMoves if applicable
    if mc_routine.clusters_freq > 0
        gc_params.cluster_p = mc_routine.initial_cluster_p
        gc_params.cluster_accepted = 0.0
        gc_params.cluster_total = 0.0
        empty!(gc_params.cluster_p_history)
        empty!(gc_params.cluster_accept_history)
        empty!(gc_params.cluster_adjust_iterations)
    end

    df = DataFrame(iter=Int[], omega=Float64[], energy=Float64[], num_particles=Int[])

    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)

        iter, omega, energy, n_par, liveset, gc_params = nested_sampling_step!(
            liveset, gc_params, mc_routine; ns_iteration=i)

        @debug "GC-NS step $i, iter: $iter, omega: $omega, energy: $energy, N: $n_par"

        if gc_params.fail_count >= gc_params.allowed_fail_count
            @warn "GC-NS: Failed $(gc_params.allowed_fail_count) times in a row."
            gc_params.fail_count = 0
        end

        if !(iter isa typeof(missing))
            push!(df, (iter, omega.val, energy.val, n_par))
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
end

"""
    IdealGasReferencedGCNSParameters(;
        mc_steps=100, reference_fugacity=1.0, energy_perturbation=1e-12,
        random_seed=1234, fail_count=0, allowed_fail_count=10,
        cluster_p=0.3, cluster_accepted=0.0, cluster_total=0.0,
        cluster_p_history=Float64[], cluster_accept_history=Float64[],
        cluster_adjust_iterations=Int[])

Convenience constructor for `IdealGasReferencedGCNSParameters`.

`reference_fugacity` (z0) must be positive. Choose z0 near the target
fugacity range `exp(βμ)` of interest: post-run reweighting to a target (μ, T)
carries a factor `(exp(βμ)/z0)^N` per sample, and its effective sample size
degrades as `|βμ − ln z0|` grows (roughly beyond `1/√Var(N)`).

The `cluster_*` fields are mutable runtime state for adaptive cluster move
tuning, initialized from the static configuration on `MCGrandCanonicalMoves`
at the start of `ideal_gas_referenced_nested_sampling` when `clusters_freq > 0`.
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
    to_walk = deepcopy(ats[parent_idx])

    # Decorrelate via GC MCMC: with mu = 0 the Ω ceiling reduces to an energy
    # ceiling, and z0 weights the insert/delete acceptance to preserve the
    # Bernoulli(z0/(1+z0)) prior
    accept, rate, to_walk, cl_accepted, cl_total = MC_grand_canonical_walk!(
        params.mc_steps, to_walk, h, emax_val, 0.0;
        p_move=mc_routine.p_move, p_insert=mc_routine.p_insert,
        energy_perturb=params.energy_perturbation,
        z0=params.reference_fugacity,
        clusters_freq=mc_routine.clusters_freq,
        swaps_freq=mc_routine.swaps_freq,
        cluster_p=params.cluster_p)

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

# Returns
- `df::DataFrame`: Columns `[:iter, :emax, :num_particles]`.
- `liveset::LatticeGasWalkers`: The final liveset (surviving walkers).
- `params::IdealGasReferencedGCNSParameters`: Updated parameters.
"""
function ideal_gas_referenced_nested_sampling(liveset::LatticeGasWalkers,
                                              params::IdealGasReferencedGCNSParameters,
                                              n_steps::Int64,
                                              mc_routine::MCGrandCanonicalMoves,
                                              save_strategy::DataSavingStrategy)
    # Initialize walkers as i.i.d. draws from the Bernoulli(z0/(1+z0)) prior
    _init_ideal_gas_ref_walkers!(liveset, params)

    # Initialize cluster_p and reset counters from MCGrandCanonicalMoves if applicable
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

        iter, emax, n_par, liveset, params = nested_sampling_step!(
            liveset, params, mc_routine; ns_iteration=i)

        @debug "IG-ref GC-NS step $i, iter: $iter, emax: $emax, N: $n_par"

        if params.fail_count >= params.allowed_fail_count
            @warn "IG-ref GC-NS: Failed $(params.allowed_fail_count) times in a row."
            params.fail_count = 0
        end

        if !(iter isa typeof(missing))
            push!(df, (iter, emax.val, n_par))
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