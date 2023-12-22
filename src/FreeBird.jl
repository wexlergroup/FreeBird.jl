module FreeBird

include("Potentials.jl")
using .Potentials
export LennardJonesParameters, lennard_jones_energy

include("AbstractWalkers.jl")
using .AbstractWalkers
export AtomsWalkers, LennardJonesAtomsWalkers

include("AtomsMCMoves.jl")
using .AtomsMCMoves
export compute_total_energy, random_walk, AtomsSystem

include("NestedSamplingLoops.jl")
using .NestedSamplingLoops
export NestedSamplingParameters
export highest_energy, sort_by_energy!, nested_sampling_step!, assign_lj_energies!

end
