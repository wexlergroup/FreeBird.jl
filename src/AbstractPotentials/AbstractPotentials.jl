"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful

export AbstractPotential
export SingleComponentParameterSets, MultiComponentParameterSets


export CompositeParameterSets

export LJParameters, lj_energy
# export CompositeLJParameters
export LennardJonesParameterSets
export GuptaParameters, gupta_attraction_squared, gupta_repulsion

abstract type AbstractPotential end

abstract type SingleComponentParameterSets <: AbstractPotential end
abstract type MultiComponentParameterSets <: AbstractPotential end

include("composite_paramsets.jl")

include("lennardjones.jl")

include("gupta.jl")

end # module Potentials