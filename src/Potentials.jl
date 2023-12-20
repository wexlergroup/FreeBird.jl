module Potentials

using Parameters

export LennardJonesParameters, lennardjonesenergy

function lennardjonesenergy(sigma::Float64, epsilon::Float64, r::Float64)
    r6 = (sigma / r)^6
    r12 = r6^2
    return 4 * epsilon * (r12 - r6)
end

struct LennardJonesParameters
    ϵ::Float64
    σ::Float64
    cutoff::Float64
    shift::Float64
end


function LennardJonesParameters(;ϵ=1.0, σ=1.0, cutoff=Inf, shift=false)
    if isa(shift,Bool) && shift
        shiftenergy = lennardjonesenergy(Float64(ϵ), Float64(σ), Float64(cutoff))
        return LennardJonesParameters(Float64(ϵ), Float64(σ), Float64(cutoff), shiftenergy)
    else
        return LennardJonesParameters(Float64(ϵ), Float64(σ), Float64(cutoff), Float64(shift))
    end
end





end # module Potentials