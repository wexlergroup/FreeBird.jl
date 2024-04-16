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