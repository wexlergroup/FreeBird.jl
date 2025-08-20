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
export GuptaParameters, gupta_attraction_squared, gupta_repulsion

abstract type AbstractPotential end

include("lennardjones.jl")

include("gupta.jl")

end # module Potentials