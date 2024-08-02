"""
Method for defining and implementing Hamiltonians.
"""
module Hamiltonians

using Unitful

export LatticeGasHamiltonian

abstract type ClassicalHamiltonian end

"""
    struct LatticeGasHamiltonian <: ClassicalHamiltonian

The `LatticeGasHamiltonian` struct represents the parameters for a lattice-gas Hamiltonian.

# Fields
- `adsorption_energy::typeof(1.0u"eV")`: The energy of adsorption on the lattice.
- `nn_interaction_energy::typeof(1.0u"eV")`: The energy of nearest-neighbor interactions.
- `nnn_interaction_energy::typeof(1.0u"eV")`: The energy of next-nearest-neighbor interactions.

"""

struct LatticeGasHamiltonian <: ClassicalHamiltonian
    adsorption_energy::typeof(1.0u"eV")
    nn_interaction_energy::typeof(1.0u"eV")
    nnn_interaction_energy::typeof(1.0u"eV")
end

function Base.show(io::IO, hamiltonian::LatticeGasHamiltonian)
    println(io, "LatticeGasHamiltonian:")
    println(io, "    adsorption_energy:      ", hamiltonian.adsorption_energy)
    println(io, "    nn_interaction_energy:  ", hamiltonian.nn_interaction_energy)
    println(io, "    nnn_interaction_energy: ", hamiltonian.nnn_interaction_energy)
end

end # module Hamiltonians