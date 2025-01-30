using FreeBird

using Statistics
using StaticArrays
using DataFrames

using Unitful
using Unitful: DimensionError
using AtomsBase
using ExtXYZ
using Test

include("test-AbstractHamiltonians.jl")

include("test-AbstractLiveSets.jl")

include("test-AbstractPotentials.jl")

include("test-AbstractWalkers.jl")

include("test-AnalysisTools.jl")

include("test-EnergyEval.jl")

include("test-MonteCarloMoves.jl")

include("test-SamplingSchemes.jl")

include("test-FreeBirdIO.jl")
