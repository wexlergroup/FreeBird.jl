module Potentials

# using Parameters
using Unitful

export LJParameters, lj_energy

"""
    struct LJParameters

A struct representing the parameters for the Lennard-Jones potential.

# Fields
- `epsilon::Float64`: The depth of the potential well. In units of eV.
- `sigma::Float64`: The distance at which the potential is zero. In units of Å.
- `cutoff::Float64`: The distance at which the potential is truncated. In units of sigma.
- `shift::Float64`: The energy shift applied to the potential calculated at the cutoff distance.
"""
struct LJParameters
    epsilon::typeof(1.0u"eV")
    sigma::typeof(1.0u"Å")
    cutoff::Float64
    shift::typeof(0.0u"eV")
end

"""
    LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
A constructor for the LJParameters struct with default values for the 
Lennard-Jones potential with no cutoff or shift. The `shift` parameter can be 
specified as a boolean, if `true`, the shift energy is calculated automatically 
at the cutoff distance; or as a `Float64`, in which case the value is used directly.

# Example

```jldoctest
julia> lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)
LJParameters(0.1 eV, 2.5 Å, 3.5, 0.0 eV)

julia> lj = LJParameters(sigma=2.5)
LJParameters(1.0 eV, 2.5 Å, Inf, 0.0 eV)

julia> lj = LJParameters(cutoff=3.5,shift=5.0)
LJParameters(1.0 eV, 1.0 Å, 3.5, 5.0 eV)

julia> lj = LJParameters(cutoff=3.5,shift=true)
LJParameters(1.0 eV, 1.0 Å, 3.5, -0.0021747803916549904 eV)

julia> lj = LJParameters(cutoff=3.5,shift=false)
LJParameters(1.0 eV, 1.0 Å, 3.5, 0.0 eV)

```
"""
function LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
    if isa(shift,Bool) && shift
        shiftenergy = lj_energy(Float64(epsilon)*u"eV", Float64(sigma)*u"Å", cutoff*Float64(sigma)*u"Å")  
    else
        shiftenergy = Float64(shift)*u"eV"
    end
    return LJParameters(epsilon*u"eV", sigma*u"Å", cutoff, shiftenergy)
end

"""
    lj_energy(epsilon::eV, sigma::Å, r::Å)

Compute the Lennard-Jones potential energy between two particles.

Arguments:
- `epsilon`: the energy scale of the potential (in electron volts)
- `sigma`: the distance scale of the potential (in angstroms)
- `r`: the distance between the particles (in angstroms)

Returns:
- The Lennard-Jones potential energy between the particles.

"""
function lj_energy(epsilon::typeof(1.0u"eV"), sigma::typeof(1.0u"Å"), r::typeof(1.0u"Å"))
    r6 = (sigma / r)^6
    r12 = r6^2
    return 4 * epsilon * (r12 - r6)
end

"""
    lj_energy(r::Å, lj::LJParameters)
Calculate the Lennard-Jones energy for a given distance `r` and LJParameters `lj`.
"""
function lj_energy(r::typeof(1.0u"Å"), lj::LJParameters)
    if r > lj.cutoff * lj.sigma
        return 0.0
    else
        return lj_energy(lj.epsilon, lj.sigma, r) - lj.shift
    end
end





end # module Potentials