"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful

export AbstractPotential
export LJParameters, lj_energy
export CompositeLJParameters
export LennardJonesParameterSets

abstract type AbstractPotential end

include("lennardjones.jl")

end # module Potentials