"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful

export AbstractPotential
export SingleComponentPotential, MultiComponentPotential

export PotentialStyle
export Pairwise, ManyBody


export CompositeParameterSets

export LJParameters, lj_energy
# export CompositeLJParameters
export LennardJonesParameterSets
export GuptaParameters #, gupta_attraction_squared, gupta_repulsion

export pair_energy, two_body_energy, many_body_energy, total_energy


abstract type AbstractPotential end

abstract type SingleComponentPotential{T} <: AbstractPotential end
abstract type MultiComponentPotential <: AbstractPotential end


abstract type PotentialStyle end

abstract type Pairwise <: PotentialStyle end
abstract type ManyBody <: PotentialStyle end

include("composite_paramsets.jl")

include("lennardjones.jl")

include("gupta.jl")

end # module Potentials