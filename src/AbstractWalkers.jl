module AbstractWalkers

# using ExtXYZ
using ..Potentials

export AtomsWalkers, LennardJonesAtomsWalkers

abstract type AtomsWalkers end

struct LennardJonesAtomsWalkers <: AtomsWalkers
    walkers::Vector
    lennard_jones_potential::LennardJonesParameters
end

end # module AbstractWalkers