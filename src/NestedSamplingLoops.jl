module NestedSamplingLoops

using ExtXYZ
using Setfield

using ..AtomsMCMoves
using ..Potentials
using ..AbstractWalkers

export NestedSamplingParameters
export highest_energy, sort_by_energy!, nested_sampling_step!, assign_lj_energies!

abstract type SamplingParameters end

struct NestedSamplingParameters <: SamplingParameters
    mc_steps::Int64
    step_size::Float64
    frozen::Int64
end

function assign_lj_energies!(liveset::LennardJonesAtomsWalkers)
    ats = liveset.walkers
    lj = liveset.lennard_jones_potential
    @debug println("ats[1].system_data.energy: ", ats[1].system_data.energy)
    # for at in ats
    #     println("at.system_data.energy: ", at.system_data.energy)
    #     at = @set at.system_data.energy = compute_total_energy(at, lj)
    # end
    for i in eachindex(ats)
        at = ats[i]
        at = @set at.system_data.energy = compute_total_energy(at, lj)
        # println("at.system_data.energy[$i]: ", at.system_data.energy)
        ats[i] = at
    end
end

function sort_by_energy!(liveset::LennardJonesAtomsWalkers)
    sort!(liveset.walkers, by = x -> x.system_data.energy, rev=true)
    # println("after sort ats[1].system_data.energy: ", ats[1].system_data.energy)
    return liveset
end



function nested_sampling_step!(liveset::LennardJonesAtomsWalkers, ns_params::NestedSamplingParameters)
    sort_by_energy!(liveset)
    @debug println("liveset.walkers[1].system_data.energy: ", liveset.walkers[1].system_data.energy)
    ats = liveset.walkers
    lj = liveset.lennard_jones_potential
    emax = liveset.walkers[1].system_data.energy::Float64
    @debug println("emax: ", emax) # debug
    to_walk = ats[1]::Atoms
    n_accept, at = random_walk(ns_params.mc_steps, to_walk, lj, ns_params.step_size, emax, ns_params.frozen)
    @debug println("n_accept: ", n_accept) # debug
    push!(ats, at)
    popfirst!(ats)
    liveset = @set liveset.walkers = ats
    return emax, liveset
end







end # module NestedSamplingLoops