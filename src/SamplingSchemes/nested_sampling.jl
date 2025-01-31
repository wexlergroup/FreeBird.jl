"""
    mutable struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in the nested sampling scheme.

# Fields
- `mc_steps::Int64`: The number of total Monte Carlo moves to perform.
- `initial_step_size::Float64`: The initial step size, which is the fallback step size if MC routine fails to accept a move.
- `step_size::Float64`: The on-the-fly step size used in the sampling process.
- `step_size_lo::Float64`: The lower bound of the step size.
- `step_size_up::Float64`: The upper bound of the step size.
- `fail_count::Int64`: The number of failed MC moves in a row.
- `allowed_fail_count::Int64`: The maximum number of failed MC moves allowed before resetting the step size.

"""
mutable struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    initial_step_size::Float64
    step_size::Float64
    step_size_lo::Float64
    step_size_up::Float64
    fail_count::Int64
    allowed_fail_count::Int64
end

mutable struct LatticeNestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    energy_perturbation::Float64
    fail_count::Int64
    allowed_fail_count::Int64
end



"""
    abstract type MCRoutine

An abstract type representing a Monte Carlo routine.
"""
abstract type MCRoutine end

"""
    struct MCRandomWalkMaxE <: MCRoutine
A type for generating a new walker by performing a random walk for decorrelation on the highest-energy walker.
"""
struct MCRandomWalkMaxE <: MCRoutine end

"""
    struct MCRandomWalkClone <: MCRoutine
A type for generating a new walker by cloning an existing walker and performing a random walk for decorrelation.
"""
struct MCRandomWalkClone <: MCRoutine end

"""
    struct MCNewSample <: MCRoutine
A type for generating a new walker from a random configuration. Currently, it is intended to use this routine for lattice gas systems.
"""
struct MCNewSample <: MCRoutine end

""" 
    struct MCMixedMoves <: MCRoutine
A type for generating a new walker by performing random walks and swapping atoms. Currently, it is intended to use this routine for
multi-component systems. The actual number of random walks and swaps to perform is determined by the weights of the fields `walks_freq` and `swaps_freq`.
For example, if `walks_freq=4` and `swaps_freq=1`, then the routine will perform 4 random walks and 1 atom swap in each cycle. The counters 
`walks_counter` and `swaps_counter` are used to keep track of the number of random walks and swaps performed in the current cycle.

# Fields
- `walks_freq::Int`: The frequency of random walks to perform.
- `swaps_freq::Int`: The frequency of atom swaps to perform.
- `walks_counter::Int`: The counter for random walks.
- `swaps_counter::Int`: The counter for atom swaps.
"""
mutable struct MCMixedMoves <: MCRoutine
    walks_freq::Int
    swaps_freq::Int
    walks_counter::Int
    swaps_counter::Int
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
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRandomWalk)

Perform a single step of the nested sampling algorithm using the Monte Carlo random walk routine.

Arguments
- `liveset::AtomWalkers`: The set of atom walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCRandomWalk`: The Monte Carlo random walk routine. Currently within this function, the random walk is only applied to the highest-energy walker, i.e., the one being culled.

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
    accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    @info "iter: $(liveset.walkers[1].iter), acceptance rate: $rate, emax: $emax, is_accepted: $accept, step_size: $(ns_params.step_size)"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        @warn "Failed to accept MC move"
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
    to_walk = deepcopy(rand(ats[2:end]))
    if mc_routine.walks_counter  + mc_routine.swaps_counter == 0
        mc_routine.walks_counter += mc_routine.walks_freq
        mc_routine.swaps_counter += mc_routine.swaps_freq
    end
    # @show mc_routine
    if mc_routine.walks_counter > 0
        accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
        mc_routine.walks_counter -= 1
        @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: $(round(ns_params.step_size; sigdigits=4))"
    elseif mc_routine.swaps_counter > 0
        accept, rate, at = MC_random_swap!(ns_params.mc_steps, to_walk, lj, emax)
        mc_routine.swaps_counter -= 1
        @info "iter: $(liveset.walkers[1].iter), acceptance rate: $(round(rate; sigdigits=4)), emax: $(round(typeof(1.0u"eV"), emax; sigdigits=10)), is_accepted: $accept, step_size: swap"
    end
    
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate)
    return iter, emax, liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers, ns_params::LatticeNestedSamplingParameters, mc_routine::MCRoutine; accept_same_config::Bool=true)

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
                               ns_params::LatticeNestedSamplingParameters, 
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

    @info "iter: $(liveset.walkers[1].iter), acceptance rate: $rate, emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    # adjust_step_size(ns_params, rate)
    return iter, emax, liveset, ns_params
end

"""
    nested_sampling_step!(liveset::LatticeGasWalkers, ns_params::LatticeNestedSamplingParameters, mc_routine::MCNewSample; accept_same_config::Bool=false)

Perform a single step of the nested sampling algorithm.

This function takes a `liveset` of lattice gas walkers, `ns_params` containing the parameters for nested sampling, and `mc_routine` representing the Monte Carlo routine for generating new samples. It performs a single step of the nested sampling algorithm by updating the liveset with a new walker.

## Arguments
- `liveset::LatticeGasWalkers`: The liveset of lattice gas walkers.
- `ns_params::LatticeNestedSamplingParameters`: The parameters for nested sampling.
- `mc_routine::MCNewSample`: The Monte Carlo routine for generating new samples.
- `accept_same_config::Bool=true`: A flag indicating whether to accept a new sample with the same configuration as an existing one.

## Returns
- `iter`: The iteration number of the liveset after the step.
- `emax`: The maximum energy of the liveset after the step.
- `liveset::LatticeGasWalkers`: The updated liveset after the step.
- `ns_params::LatticeNestedSamplingParameters`: The updated nested sampling parameters after the step.
"""
function nested_sampling_step!(liveset::LatticeGasWalkers, 
                               ns_params::LatticeNestedSamplingParameters, 
                               mc_routine::MCNewSample)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    h = liveset.hamiltonian
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,Float64} = liveset.walkers[1].energy.val

    to_walk = deepcopy(ats[1])

    accept, at = MC_new_sample!(to_walk, h, emax; energy_perturb=ns_params.energy_perturbation)

    @info "iter: $(liveset.walkers[1].iter), emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
        iter = liveset.walkers[1].iter
    else
        @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    # adjust_step_size(ns_params, rate)
    return iter, emax, liveset, ns_params
end

"""
    adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64)

Adjusts the step size of the nested sampling algorithm based on the acceptance rate. 
    If the acceptance rate is greater than 0.75, the step size is increased by 1%. 
    If the acceptance rate is less than 0.25, the step size is decreased by 1%.

# Arguments
- `ns_params::NestedSamplingParameters`: The parameters of the nested sampling algorithm.
- `rate::Float64`: The acceptance rate of the algorithm.
- `range::Tuple{Float64, Float64}`: The range of acceptance rates for adjusting the step size. Default is (0.25, 0.75).

# Returns
- `ns_params::NestedSamplingParameters`: The updated parameters with adjusted step size.
"""
function adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64; range::Tuple{Float64, Float64}=(0.25, 0.75))
    if rate > range[2] && ns_params.step_size < ns_params.step_size_up
        ns_params.step_size *= 1.05
    elseif rate < range[1] && ns_params.step_size > ns_params.step_size_lo
        ns_params.step_size *= 0.95
    end
    return ns_params
end



"""
    nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine; args...)

Perform a nested sampling loop for a given number of steps.

# Arguments
- `liveset::AtomWalkers`: The initial set of walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `n_steps::Int64`: The number of steps to perform.
- `mc_routine::MCRoutine`: The Monte Carlo routine to use.

# Returns
- `df`: A DataFrame containing the iteration number and maximum energy for each step.
- `liveset`: The updated set of walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_loop!(liveset::AtomWalkers, 
                                ns_params::NestedSamplingParameters, 
                                n_steps::Int64, 
                                mc_routine::MCRoutine,
                                save_strategy::DataSavingStrategy)
    df = DataFrame(iter=Int[], emax=Float64[])
    for i in 1:n_steps
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
        end
        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end
    return df, liveset, ns_params
end

"""
    nested_sampling_loop!(liveset::LatticeGasWalkers, ns_params::LatticeNestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine, save_strategy::DataSavingStrategy)

Perform a nested sampling loop on a lattice gas system for a given number of steps.

# Arguments
- `liveset::LatticeGasWalkers`: The initial set of walkers.
- `ns_params::LatticeNestedSamplingParameters`: The parameters for nested sampling.
- `n_steps::Int64`: The number of steps to perform.
- `mc_routine::MCRoutine`: The Monte Carlo routine to use.
- `save_strategy::DataSavingStrategy`: The strategy for saving data.

# Returns
- `df`: A DataFrame containing the iteration number and maximum energy for each step.
- `liveset`: The updated set of walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_loop!(liveset::LatticeGasWalkers,
                                ns_params::LatticeNestedSamplingParameters, 
                                n_steps::Int64, 
                                mc_routine::MCRoutine,
                                save_strategy::DataSavingStrategy)

    df = DataFrame(iter=Int[], emax=Float64[], config=Any[])
    for i in 1:n_steps
        # write_walker_every_n(liveset.walkers[1], i, save_strategy)
        
        config = liveset.walkers[1].configuration.components

        iter, emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine)
        @debug "n_step $i, iter: $iter, emax: $emax"
        if ns_params.fail_count >= ns_params.allowed_fail_count
            @warn "Failed to accept MC move $(ns_params.allowed_fail_count) times in a row. Break!"
            ns_params.fail_count = 0
            break
        end
        if !(iter isa typeof(missing))
            push!(df, (iter, emax, config))
        end
        write_df_every_n(df, i, save_strategy)
        write_ls_every_n(liveset, i, save_strategy)
    end
    return df, liveset, ns_params
end