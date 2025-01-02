using FreeBird
using Test

using StaticArrays
using DataFrames
using Unitful
using Unitful: DimensionError
using AtomsBase
using LinearAlgebra

include("test-AbstractHamiltonians.jl")

include("test-AbstractPotentials.jl")

include("test-AbstractWalkers.jl")

include("test-AnalysisTools.jl")

include("test-EnergyEval.jl")

include("test-MonteCarloMoves.jl")
