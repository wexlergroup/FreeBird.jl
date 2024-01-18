module NestedSamplingLoops

using ExtXYZ
using Setfield

using ..AtomsMCMoves
using ..Potentials
using ..AbstractWalkers
using ..EnergyEval

export NestedSamplingParameters
export sort_by_energy!, nested_sampling_step!, assign_lj_energies!
export nested_sampling_loop!

abstract type SamplingParameters end

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
    assign_lj_energies!(liveset::LJAtomWalkers)

Assigns the Lennard-Jones potential energies to the walkers in the given `liveset`.

# Arguments
- `liveset::LJAtomWalkers`: The liveset of walkers to be assigned energies.
"""
function assign_lj_energies!(liveset::LJAtomWalkers)
    ats = liveset.walkers
    lj = liveset.lj_potential
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = total_energy(at, lj)
        ats[i] = at
    end
end

function assign_lj_energies!(liveset::LJAtomWalkersWithFrozenPart)
    ats = liveset.walkers
    lj = liveset.lj_potential
    frozen = liveset.num_frozen_particles
    e_frozen = liveset.energy_frozen_particles
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = interaction_energy(at, lj; frozen=frozen) + e_frozen
        ats[i] = at
    end
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
    sort!(liveset.walkers, by = x -> x.system_data.energy, rev=true)
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
function nested_sampling_step!(liveset::AtomWalkers, ns_params::NestedSamplingParameters)
    sort_by_energy!(liveset)
    @debug println("liveset.walkers[1].system_data.energy: ", liveset.walkers[1].system_data.energy)
    ats = liveset.walkers
    lj = liveset.lj_potential
    emax = liveset.walkers[1].system_data.energy::Float64
    @debug println("emax: ", emax) # debug
    to_walk = ats[1]::Atoms
    if liveset isa LJAtomWalkersWithFrozenPart
        frozen = liveset.num_frozen_particles
        e_frozen = liveset.energy_frozen_particles
    else
        frozen = 0
        e_frozen = 0.0
    end
    n_accept, at = random_walk(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax; frozen=frozen, e_shift=e_frozen)
    @debug println("n_accept: ", n_accept) # debug
    push!(ats, at)
    popfirst!(ats)
    liveset = @set liveset.walkers = ats
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
function nested_sampling_loop!(liveset::AtomWalkers, ns_params::NestedSamplingParameters, n_steps::Int64)
    energies = zeros(Float64, n_steps)
    assign_lj_energies!(liveset)
    for i in 1:n_steps
        emax, liveset = nested_sampling_step!(liveset, ns_params)
        energies[i] = emax
    end
    return energies, liveset
end




end # module NestedSamplingLoops