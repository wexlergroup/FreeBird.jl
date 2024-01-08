module FreeBird

include("Potentials.jl")
using .Potentials
export LJParameters, lj_energy

include("AbstractWalkers.jl")
using .AbstractWalkers
export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart
export free_free_energy, frozen_energy, free_frozen_energy
export interaction_energy, total_energy

include("AtomsMCMoves.jl")
using .AtomsMCMoves
export random_walk

include("NestedSamplingLoops.jl")
using .NestedSamplingLoops
export NestedSamplingParameters
export highest_energy, sort_by_energy!, nested_sampling_step!, assign_lj_energies!
export nested_sampling_loop!

end
