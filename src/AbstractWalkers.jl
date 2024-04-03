"""
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using ..Potentials
using ..EnergyEval

export AtomWalker, AtomWalkers, LJAtomWalkers
export LatticeWalkers, Lattice2DWalker, Lattice2DWalkers
export update_walker!

abstract type AtomWalkers end

"""
    mutable struct AtomWalker

The `AtomWalker` struct represents a walker composed of atoms/molecules.

# Fields
- `configuration::FastSystem`: The configuration of the walker.
- `energy::typeof(0.0u"eV")`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.
- `num_frozen_part::Int64`: The number of frozen particles in the walker.
- `energy_frozen_part::typeof(0.0u"eV")`: The energy of the frozen particles in the walker, serves as a constant energy offset
    to the interacting part of the system.

# Constructor
```julia
AtomWalker(configuration::FastSystem; energy=0.0u"eV", iter=0, num_frozen_part=0, energy_frozen_part=0.0u"eV")
```
Create a new `AtomWalker` with the given configuration and optional energy, iteration number, number of frozen particles, and energy of the frozen particles.

"""
mutable struct AtomWalker 
    configuration::FastSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    num_frozen_part::Int64
    energy_frozen_part::typeof(0.0u"eV")
    function AtomWalker(configuration::FastSystem; energy=0.0u"eV", iter=0, num_frozen_part=0, energy_frozen_part=0.0u"eV")
        return new(configuration, energy, iter, num_frozen_part, energy_frozen_part)
    end
end

"""
    assign_lj_energies!(walkers::Vector{AtomWalker}, lj::LJParameters)

Assigns LJ energies to a list of AtomWalkers.

# Arguments
- `walkers::Vector{AtomWalker}`: A list of AtomWalkers to assign energies to.
- `lj::LJParameters`: The LJ parameters used for energy calculation.

# Returns
- `walkers::Vector{AtomWalker}`: The updated list of AtomWalkers with assigned energies.
"""
function assign_lj_energies!(walkers::Vector{AtomWalker}, lj::LJParameters)#; frozen::Int64=0, e_frozen=0.0u"eV")
    for walker in walkers
        if walker.num_frozen_part > 0
            e_frozen = frozen_energy(walker.configuration, lj, walker.num_frozen_part)
            walker.energy_frozen_part = e_frozen
        else
            e_frozen = 0.0u"eV"
        end
        e_total = interaction_energy(walker.configuration, lj; frozen=walker.num_frozen_part) + e_frozen
        walker.energy = e_total
    end
    return walkers
end

"""
    struct LJAtomWalkers <: AtomWalkers

The `LJAtomWalkers` struct contains a list of `AtomWalker` objects and the Lennard-Jones potential parameters 
    that defines the interactions between the particles in the walkers.

# Fields
- `walkers::Vector{AtomWalker}`: The list of `AtomWalkers`.
- `lj_potential::LJParameters`: The Lennard-Jones potential parameters.

# Constructors
- `LJAtomWalkers(walkers::Vector{AtomWalker}, lj_potential::LJParameters)`: Constructs a new `LJAtomWalkers` 
    object with the given walkers and Lennard-Jones potential parameters. The energies of the walkers are automatically
    assigned using the Lennard-Jones potential parameters.

"""
struct LJAtomWalkers <: AtomWalkers
    walkers::Vector{AtomWalker}
    lj_potential::LJParameters
    function LJAtomWalkers(walkers::Vector{AtomWalker}, lj_potential::LJParameters)
        assign_lj_energies!(walkers, lj_potential)
        return new(walkers, lj_potential)
    end
end

"""
    update_walker!(walker::AtomWalker, key::Symbol, value)

Update the properties of an AtomWalker object.

A convenient function that updates the value of a specific property of an AtomWalker object.

# Arguments
- `walker::AtomWalker`: The AtomWalker object to be updated.
- `key::Symbol`: The key of the property to be updated.
- `value`: The new value of the property.

# Returns
- `walker::AtomWalker`: The updated AtomWalker object.

# Example
```julia
update_walker!(walker, :energy, 10.0u"eV")
update_walker!(walker, :iter, 1)
```

"""
function update_walker!(walker::AtomWalker, key::Symbol, value)
    setproperty!(walker, key, value)
    return walker
end

"""
    LJAtomWalkers(ats::Vector{AtomWalker}, lj::LJParameters, num_frozen_part::Int64)

Create a new `LJAtomWalkers` object with the given number of frozen particles. This function is a convenient 
wrapper around the `LJAtomWalkers` constructor when the number of frozen particles is undefined or modified.

# Arguments
- `ats::Vector{AtomWalker}`: The list of `AtomWalker` objects.
- `lj::LJParameters`: The Lennard-Jones potential parameters.
- `num_frozen_part::Int64`: The number of frozen particles.

# Returns
- `LJAtomWalkers`: The `LJAtomWalkers` object with the updated number of frozen particles and energies assigned to the walkers.

"""
function LJAtomWalkers(ats::Vector{AtomWalker}, lj::LJParameters, num_frozen_part::Int64)
    for at in ats
        at.num_frozen_part = num_frozen_part
    end
    return LJAtomWalkers(ats, lj)
end

struct Lattice2DSystem
    lattice_type::Symbol  # Type of lattice (e.g., :square, :hexagonal)
    dimensions::Tuple{Int64, Int64}  # Dimensions of the lattice (rows x columns for 2D)
    site_occupancy::Matrix{Bool}  # Matrix to track whether a lattice site is occupied
    site_energy::Matrix{Float64}  # Matrix of energy values associated with each site
    external_field::Matrix{Float64}  # Matrix representing an external field affecting site energies

    # Constructor for Lattice2DSystem
    function Lattice2DSystem(lattice_type::Symbol, dims::Tuple{Int64, Int64}; external_field_strength::Float64=0.0)
        site_occupancy = falses(dims)  # Initialize all sites as unoccupied
        site_energy = zeros(dims)  # Initialize site energies to zero
        external_field = fill(external_field_strength, dims)  # Apply uniform external field
        
        # Custom initialization logic can go here, for example, setting up initial site energies based on lattice type
        
        return new(lattice_type, dims, site_occupancy, site_energy, external_field)
    end
end

abstract type LatticeWalkers end

mutable struct Lattice2DWalker
    configuration::Lattice2DSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    num_frozen_part::Int64
    energy_frozen_part::typeof(0.0u"eV")
    function Lattice2DWalker(configuration::Lattice2DSystem; energy=0.0u"eV", iter=0, num_frozen_part=0, energy_frozen_part=0.0u"eV")
        return new(configuration, energy, iter, num_frozen_part, energy_frozen_part)
    end
end

function interaction_energy(at::Lattice2DSystem, lg::LGHamiltonian; frozen::Int64=0)
    e_free_frozen = frozen == 0 ? 0.0u"eV" : free_frozen_energy(at, lg, frozen)
    return free_free_energy(at, lg; frozen=frozen) + e_free_frozen
end

function assign_lg_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
    for walker in walkers
        if walker.num_frozen_part > 0
            e_frozen = frozen_energy(walker.configuration, lg, walker.num_frozen_part)
            walker.energy_frozen_part = e_frozen
        else
            e_frozen = 0.0u"eV"
        end
        e_total = interaction_energy(walker.configuration, lg; frozen=walker.num_frozen_part) + e_frozen
        walker.energy = e_total
    end
    return walkers
end

struct Lattice2DWalkers <: LatticeWalkers
    walkers::Vector{Lattice2DWalker}
    lj_potential::LJParameters
    function Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lj_potential::LJParameters)
        assign_lj_energies!(walkers, lj_potential)
        return new(walkers, lj_potential)
    end
end

function Lattice2DWalkers(ats::Vector{Lattice2DWalker}, lj::LJParameters, num_frozen_part::Int64)
    for at in ats
        at.num_frozen_part = num_frozen_part
    end
    return Lattice2DWalkers(ats, lj)
end

end # module AbstractWalkers