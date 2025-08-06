abstract type GuptaParameterSets <: AbstractPotential end

"""
    struct GuptaParameters
The `GuptaParameters` struct represents the parameters for the Gupta potential.
#
# Fields
- `A::typeof(1.0u"eV")`: The repulsive energy scale of the potential.
- `ξ::typeof(1.0u"eV")`: The attractive energy scale of the potential.
- `p::Float64`: The exponent for the repulsive term.
- `q::Float64`: The exponent for the attractive term.
- `r0::typeof(1.0u"Å")`: The equilibrium distance for the potential.
- `cutoff::Float64`: The cutoff distance for the potential, in units of `r0`.
"""
struct GuptaParameters <: GuptaParameterSets
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

function gupta_attraction(r::typeof(1.0u"Å"), gp::GuptaParameters)
    if r > gp.cutoff * gp.r0
        return 0.0u"eV"
    else 
        return gp.ξ^2 * exp(-2 * gp.q * (r/gp.r0 - 1))
    end
end

function gupta_repulsion(r::typeof(1.0u"Å"), gp::GuptaParameters)
    if r > gp.cutoff * gp.r0
        return 0.0u"eV"
    else 
        return gp.A * exp(-gp.p * (r/gp.r0 - 1))
    end 
end