"""
    mutable struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in the nested sampling scheme.

# Fields
- `mc_steps::Int64`: The number of total Monte Carlo moves to perform.
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
            )
    NestedSamplingParameters(mc_steps, initial_step_size, step_size, step_size_lo, step_size_up, accept_range, fail_count, allowed_fail_count, energy_perturbation, random_seed)  
end

"""
    LatticeNestedSamplingParameters(;
            mc_steps::Int64=100,
            energy_perturbation::Float64=1e-12,
            fail_count::Int64=0,
            allowed_fail_count::Int64=10,
            random_seed::Int64=1234,
            )
A convenience constructor for `NestedSamplingParameters` with default values suitable for lattice systems.
"""
function LatticeNestedSamplingParameters(;
            mc_steps::Int64=100,
            energy_perturbation::Float64=1e-12,
            fail_count::Int64=0,
            allowed_fail_count::Int64=10,
            random_seed::Int64=1234,
            )
    NestedSamplingParameters(mc_steps=mc_steps, fail_count=fail_count, allowed_fail_count=allowed_fail_count, energy_perturbation=energy_perturbation, random_seed=random_seed)
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
A type for generating a new walker by performing random walks and swapping atoms. Currently, it is intended to use this routine for
multi-component systems. The actual number of random walks and swaps to perform is determined by the weights of the fields `walks_freq` and `swaps_freq`.
For example, if `walks_freq=4` and `swaps_freq=1`, then the probability of performing a random walk is 4/5, and the probability of performing a swap is 1/5.

# Fields
- `walks_freq::Int`: The frequency of random walks to perform.
- `swaps_freq::Int`: The frequency of atom swaps to perform.
"""
mutable struct MCMixedMoves <: MCRoutine
    walks_freq::Int
    swaps_freq::Int
end

"""
    struct MCRejectionSampling <: MCRoutine
A type for generating a new walker by performing rejection sampling. Currently, it is intended to use this routine for lattice gas systems.
"""
struct MCRejectionSampling <: MCRoutine end

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

Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRoutine`: The Monte Carlo routine for generating new samples. See [`MCRoutine`](@ref).

Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
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



function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutineParallel)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    if mc_routine isa MCRandomWalkMaxEParallel
        to_walk_inds = 1:nworkers()
    elseif mc_routine isa MCRandomWalkCloneParallel
        to_walk_inds = sort!(sample(2:length(ats), nworkers()))
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

function nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutineParallel)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    iter::Union{Missing,Int} = missing
    emax::Union{Vector{Missing},Vector{typeof(0.0u"eV")}} = [liveset.walkers[i].energy for i in 1:nworkers()]

    if mc_routine isa MCRandomWalkMaxEParallel
        to_walk_inds = 1:nworkers()
    elseif mc_routine isa MCRandomWalkCloneParallel
        to_walk_inds = sort!(sample(2:length(ats), nworkers()))
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

function nested_sampling_step!(liveset::LJSurfaceWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRoutine)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
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

Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCMixedMoves`: The Monte Carlo mixed moves routine.

Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCMixedMoves)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy

    # clone one of the lower energy walkers
    to_walk = deepcopy(rand(ats[2:end]))
    # determine whether to perform a random walk or a swap
    swap_prob = mc_routine.swaps_freq / (mc_routine.walks_freq + mc_routine.swaps_freq)
    
    # @show mc_routine
    if rand() > swap_prob
        accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
        @info "Swap move performed at iter: $(liveset.walkers[1].iter), accepted: $accept"
        # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: $(round(ns_params.step_size; sigdigits=4))"
    else
        accept, rate, at = MC_random_swap!(ns_params.mc_steps, to_walk, lj, emax)
        # @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: swap"
    end
    
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
                               mc_routine::MCRoutine)
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
                               mc_routine::MCNewSample)
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
                               mc_routine::MCRejectionSampling)
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
    df = DataFrame(iter=Int[], emax=Float64[])
    for i in 1:n_steps
        print_info = i % save_strategy.n_info == 0
        write_walker_every_n(liveset.walkers[1], i, save_strategy)
        iter, emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine)
        @debug "n_step $i, iter: $iter, emax: $emax"
        if ns_params.fail_count >= ns_params.allowed_fail_count
            @warn "Failed to accept MC move $(ns_params.allowed_fail_count) times in a row. Reset step size!"
            ns_params.fail_count = 0
            ns_params.step_size = ns_params.initial_step_size
        end
        if !(iter isa typeof(missing))
            push!(df, (iter, emax.val))
            if print_info
                @info "iter: $(liveset.walkers[1].iter), emax: $(emax), step_size: $(round(ns_params.step_size; sigdigits=4))"
            end
        elseif iter isa typeof(missing) && print_info
            @info "MC move failed, step: $(i), emax: $(liveset.walkers[1].energy), step_size: $(round(ns_params.step_size; sigdigits=4))"
        end
        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end
    return df, liveset, ns_params
end