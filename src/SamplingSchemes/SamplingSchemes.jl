"""
    SamplingSchemes
    
Module for defining the sampling schemes.
"""
module SamplingSchemes

using ExtXYZ
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

# Abstract types for the sampling parameters
export NestedSamplingParameters
export LatticeNestedSamplingParameters
export WangLandauParameters
export MetropolisMCParameters

# nested sampling related functions
export sort_by_energy!, nested_sampling_step!
export nested_sampling_loop!
export MCRoutine, MCRandomWalkMaxE, MCRandomWalkClone, MixedMCRoutine, MCNewSample, MCRejectionSampling

# other sampling schemes
export exact_enumeration
export wang_landau
export nvt_monte_carlo, monte_carlo_sampling


"""
    struct NestedSamplingParameters

The `NestedSamplingParameters` struct represents the parameters for various sampling algorithm.
"""
abstract type SamplingParameters end

include("nested_sampling.jl")

include("exact_enumeration.jl")

include("nvt_monte_carlo.jl")

include("wang_landau.jl")

end # module SamplingSchemes