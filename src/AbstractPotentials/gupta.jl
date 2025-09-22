abstract type GuptaParameterSets{T} <: SingleComponentPotential{T} end

"""
    struct GuptaParameters
The `GuptaParameters` struct represents the parameters for the Gupta potential (a many-body potential, also known as the second-moment approximation tight-binding potential or SMA-TB).

# Formula
The Gupta potential is given by:
```math
E_i = \\sum_j A \\exp\\left(-p \\left(\\frac{r_{ij}}{r_0} - 1\\right)\\right) - \\sqrt{\\sum_j \\xi^2 \\exp\\left(-2q \\left(\\frac{r_{ij}}{r_0} - 1\\right)\\right)}
```
where:
- ``E_i`` is the combined attraction and repulsion energy for atom ``i``.
- ``r_{ij}`` is the distance between atoms ``i`` and ``j``.
- ``A`` is the repulsive energy scale.
- ``\\xi`` is the attractive energy scale.
- ``p`` is the exponent for the repulsive term.
- ``q`` is the exponent for the attractive term.
- ``r_0`` is the nearest-neighbor distance in the bulk material.
- The potential is typically truncated at a cutoff distance, defined as a multiple of ``r_0``.
See Cleri and Rosato 1993 [Phys. Rev. B 48, 22](https://doi.org/10.1103/PhysRevB.48.22) for more details.

# Fields
- `A::typeof(1.0u"eV")`: The repulsive energy scale of the potential.
- `ξ::typeof(1.0u"eV")`: The attractive energy scale of the potential.
- `p::Float64`: The exponent for the repulsive term.
- `q::Float64`: The exponent for the attractive term.
- `r0::typeof(1.0u"Å")`: The equilibrium distance for the potential.
- `cutoff::Float64`: The cutoff distance for the potential, in units of `r0`.
"""
struct GuptaParameters <: GuptaParameterSets{ManyBody}
    A::typeof(1.0u"eV")
    ξ::typeof(1.0u"eV")
    p::Float64
    q::Float64
    r0::typeof(1.0u"Å")
    cutoff::Float64
end

"""
    GuptaParameters(;A=1.0, ξ=1.0, p=10.0, q=5.0, r0=1.0, cutoff=Inf)
A constructor for the `GuptaParameters` struct with default values for the Gupta potential with no cutoff.
"""
function GuptaParameters(;A=1.0, ξ=1.0, p=10.0, q=5.0, r0=1.0, cutoff=Inf)
    return GuptaParameters(A*u"eV", ξ*u"eV", p, q, r0*u"Å", cutoff)
end

"""
    gupta_attraction_squared(r::typeof(1.0u"Å"), gp::GuptaParameters)
Calculate the squared attractive energy term of the Gupta potential for a given interatomic distance `r`
and a set of `GuptaParameters`. It is used to define the many-body interaction in the Gupta potential.
See `many_body_energy` for more details. The square root is taken when calculating the total energy. See `total_energy`.

# Arguments
- `r::typeof(1.0u"Å")`: The interatomic distance.
- `gp::GuptaParameters`: The parameters of the Gupta potential.

# Returns
- The squared attractive energy term of the Gupta potential.
"""
function gupta_attraction_squared(r::typeof(1.0u"Å"), gp::GuptaParameters)
    if r > gp.cutoff * gp.r0
        return 0.0u"eV^2"
    else 
        return gp.ξ^2 * exp(-2 * gp.q * (r/gp.r0 - 1))
    end
end

"""
    gupta_repulsion(r::typeof(1.0u"Å"), gp::GuptaParameters)
Calculate the repulsive energy term of the Gupta potential for a given interatomic distance `r`
and a set of `GuptaParameters`. It is used to define the two-body interaction in the Gupta potential.
See `two_body_energy` for more details.

# Arguments
- `r::typeof(1.0u"Å")`: The interatomic distance.
- `gp::GuptaParameters`: The parameters of the Gupta potential.

# Returns
- The repulsive energy term of the Gupta potential.
"""
function gupta_repulsion(r::typeof(1.0u"Å"), gp::GuptaParameters)
    if r > gp.cutoff * gp.r0
        return 0.0u"eV"
    else 
        return gp.A * exp(-gp.p * (r/gp.r0 - 1))
    end 
end

"""
    two_body_energy(r::typeof(1.0u"Å"), gp::GuptaParameters)
Calculate the two-body repulsive energy between two atoms separated by a distance `r`
using the Gupta potential defined by the parameters in `gp`.
See `gupta_repulsion` for more details.
"""
two_body_energy(r::typeof(1.0u"Å"), gp::GuptaParameters) = gupta_repulsion(r, gp)

"""
    many_body_energy(r::typeof(1.0u"Å"), gp::GuptaParameters)
Calculate the squared many-body attractive energy contribution from a single atom at a distance `r`
using the Gupta potential defined by the parameters in `gp`.
See `gupta_attraction_squared` for more details.
"""
many_body_energy(r::typeof(1.0u"Å"), gp::GuptaParameters) = gupta_attraction_squared(r, gp)


"""
    total_energy(energies_two_body::Vector{<:Real}, 
                    energies_many_body::Vector{<:Real}, 
                    pot::Union{GuptaParameters,CompositeParameterSets{C, GuptaParameters}}) where C
Calculate the total energy of a system using the Gupta potential.
The total energy is computed by summing the two-body repulsive energies and subtracting the square root of the many-body attractive energies.
See `two_body_energy` and `many_body_energy` for more details.
"""
total_energy(energies_two_body::Vector{<:Real}, 
            energies_many_body::Vector{<:Real}, 
            pot::Union{GuptaParameters,CompositeParameterSets{C, GuptaParameters}}) where C = sum(energies_two_body .- sqrt.(energies_many_body)) * u"eV"