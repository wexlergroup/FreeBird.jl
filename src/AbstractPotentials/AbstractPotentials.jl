"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful
using LinearAlgebra

export AbstractPotential
export LennardJonesParametersSets
export SingleLJParametersSet
export CompositeLJParameters
export MixedParameters
export LJParameters, lj_energy
export SMD_LJParameters

abstract type AbstractPotential end

include("lennardjones.jl")

end # module Potentials