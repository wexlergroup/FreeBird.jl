"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful
using LinearAlgebra
using StaticArrays
using AtomsBase

export AbstractPotential
export PotentialParameterSets
export SingleParameterSet
export CompositeLJParameters
export MixedParameters
export LJParameters, lj_energy
export SMD_LJParameters
export pbc_dist

abstract type AbstractPotential end

include("lennardjones.jl")

end # module Potentials