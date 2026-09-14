# Changelog

## 0.4.0 — Unreleased

### Added

- Grand-canonical nested sampling for lattice systems, including ideal-gas-
  referenced and fixed-N analysis paths.
- `AtomicLattice` support across lattice walkers, moves, nested sampling and
  extXYZ serialization.
- Lazy ICET cluster-expansion energies through `ICETHamiltonian`; importing
  FreeBird still does not require the Python ICET module.
- Shared grand-canonical exact enumeration for single-component model
  lattices.
- Fixed-site lattice μVT Metropolis sampling through `MCGrandCanonicalMoves`,
  `MetropolisMCParameters` and `μvt_monte_carlo`.
- Separate C_E, C_Ω and particle-number response reporting.

### Corrected

- Nested-sampling shell weights now follow the one-based iteration convention
  and include the remaining live-set tail.
- Grand-canonical insertion and deletion use the reverse-proposal ratio needed
  for detailed balance.
- Grand-canonical samplers consume their configured random seed and terminate
  explicitly on a stalled or degenerate live set.
- Neighbor construction and periodic hollow-site matching no longer depend on
  duplicate or ASE-mutating implementations.

### Breaking changes

- `AtomicLattice` now uses an index-keyed occupation mask as its source of
  truth; its ASE object is a synchronized cache rather than the occupation
  representation.
- Grand-canonical nested-sampling rows store bare interaction energy E. Legacy
  rows containing H = E − μN must remain tagged as `grand_H_v0`; current
  analysis accepts only `bare_E_v1`.
- `MetropolisMCParameters` has an optional `chemical_potentials` field for
  lattice μVT sampling.
- The scripts-only lattice walkers, Hamiltonians, parameter types, sampling
  routines and analysis duplicates have been removed in favor of FreeBird's
  package API.

### Migration notes

- Use `MCGrandCanonicalMoves` instead of `MCGrandCanonicalLattice`.
- Use `MetropolisMCParameters(...; chemical_potentials=μ_values)` instead of
  `GCMetropolisMCParameters`.
- Use `monte_carlo_sampling(routine, lattice, h, parameters)` or the lower-
  level `μvt_monte_carlo` instead of the former scripts-only prototype.
- Lattice μVT result tables report bare mean energy in `energy`, grand-
  potential response in `c_omega`, mean coverage in `cov`, and μ in
  `chemical_potential`. The chemical-potential contribution is used in the
  acceptance rule but is not folded into the stored energy.
