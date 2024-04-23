"""
Method for defining and implementing Hamiltonians.
"""
module Hamiltonians

using Unitful

export LGHamiltonian

"""
    struct LGHamiltonian

The `LGHamiltonian` struct represents the parameters for a lattice-gas Hamiltonian.

# Fields
- `adsorption_energy::typeof(1.0u"eV")`: The energy of adsorption on the lattice.
- `nn_interaction_energy::typeof(1.0u"eV")`: The energy of nearest-neighbor interactions.
- `nnn_interaction_energy::typeof(1.0u"eV")`: The energy of next-nearest-neighbor interactions.

"""

struct LGHamiltonian
    adsorption_energy::typeof(1.0u"eV")
    nn_interaction_energy::typeof(1.0u"eV")
    nnn_interaction_energy::typeof(1.0u"eV")
end

end # module Hamiltonians