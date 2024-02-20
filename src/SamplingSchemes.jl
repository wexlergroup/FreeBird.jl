module SamplingSchemes

using ExtXYZ
using Setfield
using Unitful
using DataFrames

using ..AtomsMCMoves
using ..Potentials
using ..AbstractWalkers
using ..EnergyEval
using ..FreeBirdIO

export NestedSamplingParameters
export sort_by_energy!, nested_sampling_step!
export nested_sampling_loop!
export MCRoutine, MCRandomWalk, MCDemonWalk, MixedMCRoutine

abstract type SamplingParameters end


"""
    mutable struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in the nested sampling scheme.

# Fields
- `mc_steps::Int64`: The number of total Monte Carlo moves to perform.
- `step_size::Float64`: The step size used in the sampling process.
- `fail_count::Int64`: The number of failed MC moves in a row. Used to terminate the sampling process if it exceeds a certain threshold.

"""
mutable struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    step_size::Float64
    fail_count::Int64
end



"""
    abstract type MCRoutine

An abstract type representing a Monte Carlo routine.
"""
abstract type MCRoutine end

"""
    struct MCRandomWalk <: MCRoutine
A type representing a Monte Carlo random walk sampling scheme.
"""
struct MCRandomWalk <: MCRoutine end

"""
    MCDemonWalk <: MCRoutine

A Monte Carlo routine for performing demon walks.

# Fields
- `e_demon_tolerance::typeof(0.0u"eV")`: The tolerance for the energy difference in the demon walk.
- `demon_energy_threshold::typeof(0.0u"eV")`: The energy threshold for the demon walk.
- `demon_gain_threshold::typeof(0.0u"eV")`: The gain threshold for the demon during each walk.
- `max_add_steps::Int`: The maximum number of steps to add in the demon walk if the demon energy `e_demon` is higher than the `demon_energy_threshold`.
"""
@kwdef struct MCDemonWalk <: MCRoutine 
    e_demon_tolerance::typeof(0.0u"eV")=1e-9u"eV"
    demon_energy_threshold::typeof(0.0u"eV")=Inf*u"eV"
    demon_gain_threshold::typeof(0.0u"eV")=Inf*u"eV"
    max_add_steps::Int=1_000_000
end

"""
    struct MixedMCRoutine <: MCRoutine

A mutable struct representing a mixed Monte Carlo routine, where the main routine is used for the majority of the steps, 
    and the backup routine is used when the main routine fails to accept a move. Currently, it is intended to use `MCRandomWalk` 
    as the main routine and `MCDemonWalk` as the backup routine.

# Fields
- `main_routine::MCRoutine`: The main Monte Carlo routine.
- `back_up_routine::MCRoutine`: The backup Monte Carlo routine.
- `ns_params_main::NestedSamplingParameters`: The nested sampling parameters for the main routine.
- `ns_params_back_up::NestedSamplingParameters`: The nested sampling parameters for the backup routine.
"""
@kwdef mutable struct MixedMCRoutine <: MCRoutine
    main_routine::MCRoutine=MCRandomWalk()
    back_up_routine::MCRoutine=MCDemonWalk()
    ns_params_main::NestedSamplingParameters=NestedSamplingParameters(1000, 0.1, 0)
    ns_params_back_up::NestedSamplingParameters=NestedSamplingParameters(10000, 0.01, 0)
end




"""
    sort_by_energy!(liveset::LJAtomWalkers)

Sorts the walkers in the liveset by their energy in descending order.

# Arguments
- `liveset::LJAtomWalkers`: The liveset of walkers to be sorted.

# Returns
- `liveset::LJAtomWalkers`: The sorted liveset.
"""
function sort_by_energy!(liveset::AtomWalkers)
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
function update_iter!(liveset::AtomWalkers)
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
- `mc_routine::MCRandomWalk`: The Monte Carlo random walk routine.

Returns
- `iter`: The iteration number after the step.
- `emax`: The highest energy recorded during the step.
- `liveset`: The updated set of atom walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRandomWalk)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    to_walk = ats[1]
    accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    @info "acceptance rate: $rate, emax: $emax, is_accepted: $accept, step_size: $(ns_params.step_size)"
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
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDemonWalk)

Perform a single step of the nested sampling algorithm using the Monte Carlo demon walk routine.

# Arguments
- `liveset::AtomWalkers`: The set of atom walkers representing the current state of the system.
- `ns_params::NestedSamplingParameters`: The parameters for the nested sampling algorithm.
- `mc_routine::MCDemonWalk`: The parameters for the Monte Carlo demon walk.

# Returns
- `iter`: The iteration number after the step.
- `emax`: The maximum energy recorded during the step.
- `liveset`: The updated set of atom walkers after the step.
- `ns_params`: The updated nested sampling parameters after the step.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDemonWalk)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    iter::Union{Missing,Int} = missing
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    pick_a_random_walker = rand(ats[2:end])
    to_walk = deepcopy(pick_a_random_walker)
    accept, rate, at, _, _ = MC_nve_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size;
                                        e_demon_tolerance=mc_routine.e_demon_tolerance,
                                        demon_energy_threshold=mc_routine.demon_energy_threshold,
                                        demon_gain_threshold=mc_routine.demon_gain_threshold,
                                        max_add_steps=mc_routine.max_add_steps)
    @info "acceptance rate: $rate, emax: $emax, is_accepted: $accept"
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
    adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64)

Adjusts the step size of the nested sampling algorithm based on the acceptance rate. 
    If the acceptance rate is greater than 0.75, the step size is increased by 1%. 
    If the acceptance rate is less than 0.25, the step size is decreased by 1%.

# Arguments
- `ns_params::NestedSamplingParameters`: The parameters of the nested sampling algorithm.
- `rate::Float64`: The acceptance rate of the algorithm.

# Returns
- `ns_params::NestedSamplingParameters`: The updated parameters with adjusted step size.
"""
function adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64)
    if rate > 0.75
        ns_params.step_size *= 1.01
    elseif rate < 0.25
        ns_params.step_size *= 0.99
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

# Keyword Arguments
- `args...`: Additional arguments.

# Returns
- `df`: A DataFrame containing the iteration number and maximum energy for each step.
- `liveset`: The updated set of walkers.
- `ns_params`: The updated nested sampling parameters.
"""
function nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine; args...)
    df = DataFrame(iter=Int[], emax=typeof(0.0u"eV")[])
    for i in 1:n_steps
        iter, emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine)
        if ns_params.fail_count >= 10
            @warn "Failed to accept MC move 10 times in a row. Break!"
            ns_params.fail_count = 0
            break
        end
        if !(iter isa typeof(missing))
            push!(df, (iter, emax))
        end
    end
    return df, liveset, ns_params
end


"""
    nested_sampling_loop!(liveset::AtomWalkers, n_steps::Int64, mc_routine::MixedMCRoutine, save_strategy::DataSavingStrategy)

Perform a nested sampling loop for a given number of steps.

# Arguments
- `liveset::AtomWalkers`: The initial set of walkers.
- `n_steps::Int64`: The number of steps to perform.
- `mc_routine::MixedMCRoutine`: The mixed Monte Carlo routine to use.
- `save_strategy::DataSavingStrategy`: The strategy for saving data.

# Returns
- `df::DataFrame`: The data frame containing the iteration number and maximum energy for each step.
- `liveset::AtomWalkers`: The updated set of walkers.
- `mc_routine.ns_params_main`: The updated nested sampling parameters for the main routine.
"""
function nested_sampling_loop!(liveset::AtomWalkers,  
                                n_steps::Int64, 
                                mc_routine::MixedMCRoutine,
                                save_strategy::DataSavingStrategy)
    df = DataFrame(iter=Int[], emax=typeof(0.0u"eV")[])
    for i in 1:n_steps
        iter, emax, liveset, mc_routine.ns_params_main = nested_sampling_step!(liveset, mc_routine.ns_params_main, mc_routine.main_routine)
        if mc_routine.ns_params_main.fail_count >= 10
            @warn "Failed to accept $(mc_routine.main_routine) move 10 times in a row. Switching to back up routine $(mc_routine.back_up_routine)!"
            mc_routine.ns_params_main.fail_count = 0
            iter, emax, liveset, mc_routine.ns_params_back_up = nested_sampling_step!(liveset, mc_routine.ns_params_back_up, mc_routine.back_up_routine)
        end
        if !(iter isa typeof(missing))
            push!(df, (iter, emax))
        end
        write_df_every_n(df, i, save_strategy)
    end
    return df, liveset, mc_routine.ns_params_main
end





end # module NestedSamplingLoops