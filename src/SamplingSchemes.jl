module SamplingSchemes

using ExtXYZ
using Setfield
using Unitful

using ..AtomsMCMoves
using ..Potentials
using ..AbstractWalkers
using ..EnergyEval

export NestedSamplingParameters
export sort_by_energy!, nested_sampling_step!
export nested_sampling_loop!
export MCRoutine, MCRandomWalk, MCDemonWalk, MixedMCRoutine

abstract type SamplingParameters end

"""
    struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in nested sampling.

# Fields
- `mc_steps::Int64`: The number of Monte Carlo steps.
- `step_size::Float64`: The step size for each Monte Carlo step.
"""
mutable struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    step_size::Float64
    fail_count::Int64
end



abstract type MCRoutine end

struct MCRandomWalk <: MCRoutine end

@kwdef struct MCDemonWalk <: MCRoutine 
    e_demon_tolerance::typeof(0.0u"eV")=1e-9u"eV"
    demon_energy_threshold::typeof(0.0u"eV")=Inf*u"eV"
    demon_gain_threshold::typeof(0.0u"eV")=Inf*u"eV"
    max_add_steps::Int=1_000_000
end

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

function update_iter!(live_set::AtomWalkers)
    for at in live_set.walkers
        at.iter += 1
    end
end


"""
    nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters)

Performs a single step of the nested sampling algorithm.

# Arguments
- `liveset::LJAtomWalkers`: The current set of walkers in the nested sampling algorithm.
- `ns_params::NestedSamplingParameters`: The parameters for the nested sampling algorithm.

# Returns
- `emax`: The maximum energy among the walkers in the liveset.
- `liveset`: The updated liveset after performing the nested sampling step.
"""
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCRandomWalk)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    to_walk = ats[1]
    accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    @info "acceptance rate: $rate, emax: $emax, is_accepted: $accept, step_size: $(ns_params.step_size)"
    if accept
        push!(ats, at)
        popfirst!(ats)
        update_iter!(liveset)
        ns_params.fail_count = 0
    else
        @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate)
    return emax, liveset, ns_params
end

function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDemonWalk)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
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
    else
        @warn "Failed to accept MC move"
        emax = missing
        ns_params.fail_count += 1
    end
    adjust_step_size(ns_params, rate)
    return emax, liveset, ns_params
end

function adjust_step_size(ns_params::NestedSamplingParameters, rate::Float64)
    if rate > 0.9
        ns_params.step_size *= 1.05
    elseif rate < 0.1
        ns_params.step_size *= 0.95
    end
    return ns_params
end


"""
    nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64)

Perform nested sampling loop for a given number of steps.

# Arguments
- `liveset::LJAtomWalkers`: The initial set of walkers.
- `ns_params::NestedSamplingParameters`: The parameters for nested sampling.
- `n_steps::Int64`: The number of steps to perform.

# Returns
- `energies`: An array of energies at each step.
- `liveset`: The final set of walkers.
"""
function nested_sampling_loop!(liveset::AtomWalkers, 
                                ns_params::NestedSamplingParameters, 
                                n_steps::Int64, 
                                mc_routine::MCRoutine;
                                args...)
    energies = Array{Union{typeof(0.0u"eV"),Missing}}(undef, n_steps)
    for i in 1:n_steps
        emax, liveset, ns_params = nested_sampling_step!(liveset, ns_params, mc_routine)
        if ns_params.fail_count >= 10
            @warn "Failed to accept MC move 10 times in a row. Break!"
            ns_params.fail_count = 0
            break
        end
        energies[i] = emax
    end
    return energies, liveset, ns_params
end


function nested_sampling_loop!(liveset::AtomWalkers,  
                                n_steps::Int64, 
                                mc_routine::MixedMCRoutine)
    energies = Array{Union{typeof(0.0u"eV"),Missing}}(undef, n_steps)
    for i in 1:n_steps
        emax, liveset, mc_routine.ns_params_main = nested_sampling_step!(liveset, mc_routine.ns_params_main, mc_routine.main_routine)
        if mc_routine.ns_params_main.fail_count >= 10
            @warn "Failed to accept $(mc_routine.main_routine) move 10 times in a row. Switching to back up routine $(mc_routine.back_up_routine)!"
            mc_routine.ns_params_main.fail_count = 0
            emax, liveset, mc_routine.ns_params_back_up = nested_sampling_step!(liveset, mc_routine.ns_params_back_up, mc_routine.back_up_routine)
        end
        energies[i] = emax
    end
    return energies, liveset, mc_routine.ns_params_main
end





end # module NestedSamplingLoops