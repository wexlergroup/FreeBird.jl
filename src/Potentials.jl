module Potentials

using Parameters

export LennardJonesParameters, lennard_jones_energy

"""
    struct LennardJonesParameters

A struct representing the parameters for the Lennard-Jones potential.

# Fields
- `epsilon::Float64`: The depth of the potential well.
- `sigma::Float64`: The distance at which the potential is zero.
- `cutoff::Float64`: The distance at which the potential is truncated.
- `shift::Float64`: The energy shift applied to the potential caclulated at the cutoff distance.
"""
struct LennardJonesParameters
    epsilon::Float64
    sigma::Float64
    cutoff::Float64
    shift::Float64
end

"""
    LennardJonesParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=false)
A constructor for the LennardJonesParameters struct with default values for the 
Lennard-Jones potential with no cutoff or shift. The `shift` parameter can be 
specified as a boolean, if `true`, the shift energy is calculated automatically 
at the cutoff distance; or as a `Float64`, in which case the value is used directly.

# Example

```jldoctest
julia> lj = LennardJonesParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)
LennardJonesParameters(0.1, 2.5, 3.5, 0.0)

julia> lj = LennardJonesParameters(sigma=2.5)
LennardJonesParameters(1.0, 2.5, Inf, 0.0)

julia> lj = LennardJonesParameters(cutoff=3.5,shift=5.0)
LennardJonesParameters(1.0, 1.0, 3.5, 5.0)

julia> lj = LennardJonesParameters(cutoff=3.5,shift=true)
LennardJonesParameters(1.0, 1.0, 3.5, -0.0021747803916549904)

julia> lj = LennardJonesParameters(cutoff=3.5,shift=false)
LennardJonesParameters(1.0, 1.0, 3.5, 0.0)

```
"""
function LennardJonesParameters(;epsilon=1.0, sigma=1.0, cutoff=Inf, shift=false)
    if isa(shift,Bool) && shift
        shiftenergy = lennard_jones_energy(Float64(epsilon), Float64(sigma), Float64(cutoff))
        
    else
        shiftenergy = Float64(shift)
    end
    return LennardJonesParameters(Float64(epsilon), Float64(sigma), Float64(cutoff), shiftenergy)
end

function lennard_jones_energy(sigma::Float64, epsilon::Float64, r::Float64)
    r6 = (sigma / r)^6
    r12 = r6^2
    return 4 * epsilon * (r12 - r6)
end

"""
    lennard_jones_energy(r::Float64, lj::LennardJonesParameters)
Calculate the Lennard-Jones energy for a given distance `r` and LennardJonesParameters `lj`.
"""
function lennard_jones_energy(r::Float64, lj::LennardJonesParameters)
    if r > lj.cutoff
        return 0.0
    else
        return lennard_jones_energy(lj.sigma, lj.epsilon, r) - lj.shift
    end
end





end # module Potentials