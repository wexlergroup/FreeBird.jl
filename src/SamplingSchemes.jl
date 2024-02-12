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
export MCRoutine, MCRandomWalk, MCDemonWalk

abstract type SamplingParameters end

abstract type MCRoutine end

struct MCRandomWalk <: MCRoutine end

struct MCDemonWalk <: MCRoutine end

"""
    struct NestedSamplingParameters <: SamplingParameters

The `NestedSamplingParameters` struct represents the parameters used in nested sampling.

# Fields
- `mc_steps::Int64`: The number of Monte Carlo steps.
- `step_size::Float64`: The step size for each Monte Carlo step.
"""
struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    step_size::Float64
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
    emax::typeof(0.0u"eV") = liveset.walkers[1].energy
    to_walk = ats[1]
    accept, rate, at = MC_random_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax)
    @info "acceptance rate: $rate, emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
    end
    return emax, liveset
end

function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, mc_routine::MCDemonWalk)
    sort_by_energy!(liveset)
    ats = liveset.walkers
    lj = liveset.lj_potential
    emax::Union{Missing,typeof(0.0u"eV")} = liveset.walkers[1].energy
    pick_a_random_walker = rand(ats[2:end])
    to_walk = deepcopy(pick_a_random_walker)
    accept, rate, at, _, _ = MC_nve_walk!(ns_params.mc_steps, to_walk, lj, ns_params.step_size)
    @info "acceptance rate: $rate, emax: $emax, is_accepted: $accept"
    if accept
        push!(ats, at)
        popfirst!(ats)
    else
        @warn "Failed to accept MC move"
        emax = missing
    end
    return emax, liveset
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
function nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64, mc_routine::MCRoutine)
    energies = Array{Union{typeof(0.0u"eV"),Missing}}(undef, n_steps)
    for i in 1:n_steps
        emax, liveset = nested_sampling_step!(liveset, ns_params, mc_routine)
        energies[i] = emax
    end
    return energies, liveset
end




end # module NestedSamplingLoops