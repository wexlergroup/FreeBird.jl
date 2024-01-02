module AbstractWalkers

# using ExtXYZ
using ..Potentials

export AtomsWalkers, LennardJonesAtomsWalkers

abstract type AtomsWalkers end

struct LennardJonesAtomsWalkers{T} <: AtomsWalkers
    walkers::Vector{T}
    lennard_jones_potential::LennardJonesParameters
end

end # module AbstractWalkers