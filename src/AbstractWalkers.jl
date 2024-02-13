module AbstractWalkers

using AtomsBase
using Unitful
using ..Potentials
using ..EnergyEval

export AtomWalker, AtomWalkers, LJAtomWalkers
export update_walker!

abstract type AtomWalkers end

mutable struct AtomWalker 
    configuration::FastSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    num_frozen_part::Int64
    energy_frozen_part::typeof(0.0u"eV")
    function AtomWalker(configuration::FastSystem; energy=0.0u"eV", iter=0, num_frozen_part=0, energy_frozen_part=0.0u"eV")
        return new(configuration, energy, iter, num_frozen_part, energy_frozen_part)
    end
end

function assign_lj_energies!(walkers::Vector{AtomWalker}, lj::LJParameters)#; frozen::Int64=0, e_frozen=0.0u"eV")
    for walker in walkers
        if walker.num_frozen_part > 0
            e_frozen = frozen_energy(walker.configuration, lj, walker.num_frozen_part)
            walker.energy_frozen_part = e_frozen
        else
            e_frozen = 0.0u"eV"
        end
        e_total = interaction_energy(walker.configuration, lj; frozen=walker.num_frozen_part) + e_frozen
        walker.energy = e_total
    end
    return walkers
end

struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker}
    lj_potential::LJParameters
    function LJAtomWalkers(walkers::Vector{AtomWalker}, lj_potential::LJParameters)
        assign_lj_energies!(walkers, lj_potential)
        return new(walkers, lj_potential)
    end
end

function update_walker!(walker::AtomWalker, key::Symbol, value)
    setproperty!(walker, key, value)
    return walker
end

end # module AbstractWalkers