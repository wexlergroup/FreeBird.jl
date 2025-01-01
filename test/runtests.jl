using FreeBird
using Test

using StaticArrays
using DataFrames
using Unitful
using AtomsBase
using LinearAlgebra

include("test-AbstractHamiltonians.jl")

include("test-AbstractWalkers.jl")

include("test-AnalysisTools.jl")

include("test-EnergyEval.jl")

include("test-Potentials.jl")
