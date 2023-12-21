module NestedSamplingLoops

using ExtXYZ
using Setfield

using ..AtomsMCMoves
using ..Potentials
using ..AbstractWalkers

export highest_energy, sort_by_energy!, nested_sampling_step!, assign_lj_energies!

function assign_lj_energies!(liveset::LennardJonesAtomsWalkers)
    ats = liveset.walkers
    lj = liveset.lennard_jones_potential
    println("ats[1].system_data.energy: ", ats[1].system_data.energy)
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



function nested_sampling_step!(liveset::LennardJonesAtomsWalkers,mc_steps::Int64,step_size::Float64,frozen::Int64=0)
    sort_by_energy!(liveset)
    println("liveset.walkers[1].system_data.energy: ", liveset.walkers[1].system_data.energy)
    ats = liveset.walkers
    lj = liveset.lennard_jones_potential
    emax = liveset.walkers[1].system_data.energy::Float64
    println("emax: ", emax)
    to_walk = ats[1]::Atoms
    n_accept, at = random_walk(mc_steps, to_walk, lj, step_size, emax, frozen)
    println("n_accept: ", n_accept)
    push!(ats, at)
    popfirst!(ats)
    liveset = @set liveset.walkers = ats
    return emax, liveset
end



end # module NestedSamplingLoops