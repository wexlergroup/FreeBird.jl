# Custom Hamiltonians

The lattice sampling stack accepts user-defined Hamiltonians through a small extension contract: `interacting_energy(lattice, hamiltonian)` is the only place the stack evaluates a Hamiltonian, so one method serves `LatticeGasWalkers` construction, all of the nested-sampling loops, and `exact_enumeration`, on every single-component lattice geometry. Externally defined energies on fixed site sets are standard practice, machine-learned potentials among them, and this page records what a custom type must provide.

## The contract

1. Subtype `ClassicalHamiltonian`. Fields are yours, including caches: nothing in the stack copies or introspects the Hamiltonian.
2. Define one energy method returning an energy-dimensioned `Unitful.Quantity` (eV-convertible; a molar energy such as kJ/mol carries a different dimension and is rejected). A method returning a plain number raises a descriptive `ArgumentError` at walker setup. The validation lives on the walker-setup surface; the raw-lattice `wang_landau` and `nvt_monte_carlo` entry points evaluate the Hamiltonian directly and sit outside it.
3. Optionally, opt into two integrations: `AbstractHamiltonians._coupling_type(::MyHam)` to compose under `SiteFieldLatticeHamiltonian` (its absence produces a self-explanatory `ArgumentError`), and `AbstractLiveSets._n_coupled_shells(::MyHam)` to participate in the uncoupled-shell warning at setup.

A worked example, a coordination-dependent energy no shipped pair, cluster, or site-field type expresses directly:

````julia
using FreeBird
using Unitful

struct CoordinationHamiltonian <: ClassicalHamiltonian
    j::Float64
end

import FreeBird.EnergyEval: interacting_energy

function interacting_energy(lattice::SLattice, h::CoordinationHamiltonian)
    occ = lattice.components[1]
    doubly = count(s -> occ[s] &&
                        count(n -> occ[n], lattice.neighbors[s][1]) == 2,
                   eachindex(occ))
    return h.j * doubly^2 * u"eV"
end
````

With that in place, the type is a drop-in for the shipped Hamiltonians:

````julia
walkers = [LatticeWalker(deepcopy(lattice), energy=0.0u"eV", iter=0) for _ in 1:32]
ls = LatticeGasWalkers(walkers, CoordinationHamiltonian(0.1))
df, ls, params = grand_canonical_nested_sampling(ls, gc_params, 1000, mc_routine, save_strategy)
````

and `exact_enumeration(lattice, CoordinationHamiltonian(0.1))` enumerates its fixed-particle-number spectrum.

## Caching expensive energies

The evaluation path is serial (walker setup is a plain loop, and the nested-sampling loops evaluate one proposal at a time), so a plain `Dict` memoization keyed by the occupation vector is safe:

````julia
struct CachedHamiltonian <: ClassicalHamiltonian
    cache::Dict{Vector{Bool},typeof(1.0u"eV")}
end

function interacting_energy(lattice::SLattice, h::CachedHamiltonian)
    get!(h.cache, copy(lattice.components[1])) do
        expensive_energy(lattice)
    end
end
````

Two practical notes. First, define `Base.show` for cache-carrying types: the default struct `show` prints the entire cache whenever a liveset is displayed. Second, the serial-evaluation guarantee is what makes the unsynchronized `Dict` safe; wrap the cache before using any future parallel evaluation path.
