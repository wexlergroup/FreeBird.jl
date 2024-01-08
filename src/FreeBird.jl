module FreeBird

include("Potentials.jl")
using .Potentials
export LJParameters, lj_energy

include("AbstractWalkers.jl")
using .AbstractWalkers
export AtomWalkers, LJAtomWalkers, LJAtomWalkersWithFrozenPart

include("AtomsMCMoves.jl")
using .AtomsMCMoves
export compute_total_energy, random_walk

include("NestedSamplingLoops.jl")
using .NestedSamplingLoops
export NestedSamplingParameters
export highest_energy, sort_by_energy!, nested_sampling_step!, assign_lj_energies!
export nested_sampling_loop!

end
