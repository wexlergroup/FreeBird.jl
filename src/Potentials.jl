module Potentials

using Parameters

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
    epsilon::Float64
    sigma::Float64
    cutoff::Float64
    shift::Float64
end

"""
    LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=false)
A constructor for the LJParameters struct with default values for the 
Lennard-Jones potential with no cutoff or shift. The `shift` parameter can be 
specified as a boolean, if `true`, the shift energy is calculated automatically 
at the cutoff distance; or as a `Float64`, in which case the value is used directly.

# Example

```jldoctest
julia> lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)
LJParameters(0.1, 2.5, 3.5, 0.0)

julia> lj = LJParameters(sigma=2.5)
LJParameters(1.0, 2.5, Inf, 0.0)

julia> lj = LJParameters(cutoff=3.5,shift=5.0)
LJParameters(1.0, 1.0, 3.5, 5.0)

julia> lj = LJParameters(cutoff=3.5,shift=true)
LJParameters(1.0, 1.0, 3.5, -0.0021747803916549904)

julia> lj = LJParameters(cutoff=3.5,shift=false)
LJParameters(1.0, 1.0, 3.5, 0.0)

```
"""
function LJParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=true)
    if isa(shift,Bool) && shift
        shiftenergy = lj_energy(Float64(epsilon), Float64(sigma), Float64(cutoff))
        
    else
        shiftenergy = Float64(shift)
    end
    return LJParameters(Float64(epsilon), Float64(sigma), Float64(cutoff), shiftenergy)
end

function lj_energy(sigma::Float64, epsilon::Float64, r::Float64)
    r6 = (sigma / r)^6
    r12 = r6^2
    return 4 * epsilon * (r12 - r6)
end

"""
    lj_energy(r::Float64, lj::LJParameters)
Calculate the Lennard-Jones energy for a given distance `r` and LJParameters `lj`.
"""
function lj_energy(r::Float64, lj::LJParameters)
    if r > lj.cutoff * lj.sigma
        return 0.0
    else
        return lj_energy(lj.sigma, lj.epsilon, r) - lj.shift
    end
end





end # module Potentials