"""
    AbstractWalkers
    
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using LinearAlgebra
using Statistics

# AtomicLattice is backed by an ASE Atoms object. `Py` is a *field type* on the
# struct, so it is resolved when the module is defined, not when a method runs:
# without this import the package does not load at all, it fails with
# `UndefVarError: Py not defined in FreeBird.AbstractWalkers`. The other uses
# (`pyconvert` in get_positions and get_adsorbate_indicies, `ase.build.fcc100`
# and `ase.build.add_adsorbate` in the constructor and add_adsorbates!,
# `ase.visualize.view` in shows.jl) sit inside function bodies and would instead
# have failed at call time.
#
# `..AbstractPotentials` imports both of these but does not re-export `Py`, so
# it cannot supply them. Same pairing as AbstractPotentials.jl:31-32.
using ASEconvert
using PythonCall

using ..AbstractPotentials
using ..AbstractHamiltonians

export AbstractWalker
export AbstractLattice
export AtomWalker
export sort_components_by_atomic_number
export split_components
export split_components_by_chemical_species
export check_num_components
export LatticeWalker
export LatticeGeometry, SquareLattice, TriangularLattice, GenericLattice
export MLattice, SLattice, GLattice, AtomicLattice
export update_walker!
export num_sites, occupied_site_count
export view_structure

abstract type AbstractWalker end

include("atomistic_walkers.jl")

include("lattice_walkers.jl")

include("helpers.jl")

include("shows.jl")

end # module AbstractWalkers