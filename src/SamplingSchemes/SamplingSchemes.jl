"""
    SamplingSchemes
    
Module for defining the sampling schemes.
"""
module SamplingSchemes

using ExtXYZ
using Setfield
using Unitful
using DataFrames
using Combinatorics
using Random
using Statistics

using ..MonteCarloMoves
using ..AbstractPotentials
using ..AbstractWalkers
using ..AbstractLiveSets
using ..EnergyEval
using ..FreeBirdIO
using ..AbstractHamiltonians

export NestedSamplingParameters, LatticeNestedSamplingParameters
export sort_by_energy!, nested_sampling_step!
export nested_sampling_loop!
export MCRoutine, MCRandomWalkMaxE, MCRandomWalkClone, MCDemonWalk, MixedMCRoutine, MCNewSample

export exact_enumeration, nvt_monte_carlo, wang_landau

abstract type SamplingParameters end

include("nested_sampling.jl")

include("exact_enumeration.jl")

include("nvt_monte_carlo.jl")

include("wang_landau.jl")

end # module SamplingSchemes