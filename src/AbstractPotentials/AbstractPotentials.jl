"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful
using LinearAlgebra

export AbstractPotential
export LJParameters, lj_energy
export CompositeLJParameters
export LennardJonesParametersSets
export SingleLJParametersSet
export SMD_LJParameters

abstract type AbstractPotential end

include("lennardjones.jl")

end # module Potentials