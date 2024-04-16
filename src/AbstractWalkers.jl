"""
Module containing abstract definitions for walkers.
"""
module AbstractWalkers

using AtomsBase
using Unitful
using Random
using Combinatorics
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

# Lattice gas

struct Lattice2DSystem
    lattice_type::Symbol  # Type of lattice (e.g., :square, :hexagonal)
    dimensions::Tuple{Int64, Int64}  # Dimensions of the lattice (rows x columns for 2D)
    num_occ_sites::Int64  # Number of occupied sites
    site_occupancy::Matrix{Bool}  # Matrix to track whether a lattice site is occupied

    # Constructor to initialize with specified site occupancy
    function Lattice2DSystem(lattice_type::Symbol, site_occupancy::Matrix{Bool})
        dims = size(site_occupancy)
        num_occ_sites = sum(site_occupancy)
        return new(lattice_type, dims, num_occ_sites, site_occupancy)
    end

    # Constructor to initialize with dimensions and random site occupancy
    function Lattice2DSystem(lattice_type::Symbol, dims::Tuple{Int64, Int64}, num_occ_sites::Int64)
        total_sites = prod(dims)  # Total number of sites
        occupancy = vcat(fill(true, num_occ_sites), fill(false, total_sites - num_occ_sites))  # Initialize occupancy array
        shuffle!(occupancy)  # Shuffle to randomize occupied sites
        site_occupancy = reshape(occupancy, dims)  # Reshape back to matrix form
        return new(lattice_type, dims, num_occ_sites, site_occupancy)
    end
end

abstract type LatticeWalkers end

mutable struct Lattice2DWalker
    configuration::Lattice2DSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    function Lattice2DWalker(configuration::Lattice2DSystem; energy=0.0u"eV", iter=0)
        return new(configuration, energy, iter)
    end
end

struct LGHamiltonian
    adsorption_energy::typeof(1.0u"eV")
    nn_interaction_energy::typeof(1.0u"eV")
    nnn_interaction_energy::typeof(1.0u"eV")
end

function neighbors(lattice_type::Symbol, dims::Tuple{Int64, Int64}, i::Int64, j::Int64)
    if lattice_type == :square
        return [(mod(i + dx - 1, dims[1]) + 1, mod(j + dy - 1, dims[2]) + 1) 
                for (dx, dy) in [(0, 1), (0, -1), (1, 0), (-1, 0)]]
    elseif lattice_type == :hexagonal
        even_row_neighbors = [(0, -1), (0, 1), (-1, -1), (-1, 0), (1, -1), (1, 0)]
        odd_row_neighbors = [(0, -1), (0, 1), (-1, 0), (-1, 1), (1, 0), (1, 1)]
        offsets = iseven(i) ? even_row_neighbors : odd_row_neighbors
        return [(mod(i + dx - 1, dims[1]) + 1, mod(j + dy - 1, dims[2]) + 1) for (dx, dy) in offsets]
    else
        # TODO: Implement other lattice types
        error("Invalid lattice type")
    end
end

function interaction_energy(at::Lattice2DSystem, lg::LGHamiltonian)  # TODO: Generalize to arbitrary lattice
    # Nearest-neighbor interaction energy
    e_nn = 0.0u"eV"

    # If nearest-neighbor interaction energy is zero, return immediately
    if lg.nn_interaction_energy == 0.0u"eV"
        return e_nn
    end

    # Compute nearest-neighbor interaction energy
    for i in 1:at.dimensions[1]
        for j in 1:at.dimensions[2]
            if at.site_occupancy[i, j]
                for (i2, j2) in neighbors(at.lattice_type, at.dimensions, i, j)
                    e_nn += at.site_occupancy[i2, j2] ? lg.nn_interaction_energy / 2 : 0.0u"eV"
                end
            end
        end
    end

    # Next-nearest-neighbor interaction energy
    e_nnn_interaction = 0.0u"eV"

    # If next-nearest-neighbor interaction energy is zero, return immediately
    if lg.nnn_interaction_energy == 0.0u"eV"
        return e_nn
    elseif at.lattice_type == :hexagonal
        # TODO: Implement hexagonal lattice next-nearest-neighbor interaction
        error("Hexagonal lattice next-nearest-neighbor interaction not implemented yet")
    end

    # Compute next-nearest-neighbor interaction energy
    for i in 1:at.dimensions[1]
        for j in 1:at.dimensions[2]
            if at.site_occupancy[i, j]
                for (dx, dy) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]
                    i2 = i + dx
                    j2 = j + dy
                    if 1 <= i2 <= at.dimensions[1] && 1 <= j2 <= at.dimensions[2]
                        if at.site_occupancy[i2, j2]
                            e_nnn_interaction += lg.nnn_interaction_energy / 2
                        end
                    end
                end
            end
        end
    end
    return e_nn + e_nnn_interaction
end

function assign_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
    for walker in walkers
        e_adsorption = walker.configuration.num_occ_sites * lg.adsorption_energy
        e_interaction = interaction_energy(walker.configuration, lg)  # nearest-neighbor and next-nearest-neighbor interactions
        e_total = e_adsorption + e_interaction
        walker.energy = e_total
    end
    return walkers
end

struct Lattice2DWalkers <: LatticeWalkers
    walkers::Vector{Lattice2DWalker}
    lg::LGHamiltonian
    function Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
        assign_energies!(walkers, lg)
        return new(walkers, lg)
    end
end

# Lattice gas tests

function test_lattice2d_system()
    # Create four 4x4 square lattices with 10 occupied sites
    lattices = [Lattice2DSystem(:square, (4, 4), 10) for i in 1:4] 

    # Create a Lattice2DWalkers object with the lattices and a Hamiltonian
    lg = LGHamiltonian(1.0u"eV", 1.0u"eV", 1.0u"eV")
    walkers = [Lattice2DWalker(lattice) for lattice in lattices]
    lg_walkers = Lattice2DWalkers(walkers, lg)
    return lg_walkers
end

# Enumerate all possible configurations of an LxL square lattice with N occupied sites

function enumerate_lattice_configs(L::Int64, N::Int64, lattice_type::Symbol=:square)  # TODO: Throw error if N is odd for hexagonal lattice
    # Generate a list of all grid points
    grid_points = [(i, j) for i in 1:L for j in 1:L]

    # Generate all combinations of N points from the grid
    all_configs = collect(combinations(grid_points, N))

    # Convert each configuration to a site occupancy matrix
    all_configs = [reshape([in((i, j), config) for i in 1:L, j in 1:L], L, L) for config in all_configs]
    
    # Generate Lattice2DSystem objects for each configuration
    lattices = [Lattice2DSystem(lattice_type, config) for config in all_configs]

    # Compute the energy of each configuration
    if lattice_type == :square
        lg = LGHamiltonian(-0.04136319965u"eV", -0.01034079991u"eV", -(15 / 64) * 0.01034079991u"eV")
    else  # Hexagonal lattice
        lg = LGHamiltonian(-0.03102239974u"eV", -0.01034079991u"eV", 0.0u"eV")
    end
    walkers = [Lattice2DWalker(lattice) for lattice in lattices]
    lg_walkers = Lattice2DWalkers(walkers, lg)

    # Generate a list of energies
    energies = [walker.energy for walker in lg_walkers.walkers]

    return energies
end

function compute_internal_energy_versus_temperature(L::Int64, N::Int64, T_min::typeof(1.0u"K"), T_max::typeof(100.0u"K"), num_points::Int64, lattice_type::Symbol=:square)
    # Generate all possible configurations with 8 occupied sites
    energies = enumerate_lattice_configs(L, N, lattice_type)

    # Compute energy relative to the lowest energy
    minimum_energy = minimum(energies)
    energies = energies .- minimum_energy

    # Compute the partition function
    BoltzmannConstant = 8.617_333_262e-5u"eV/K"
    temperatures = range(T_min, T_max, length=num_points)
    internal_energies = [sum(energy * exp(-energy / (BoltzmannConstant * T)) for energy in energies) / sum(exp(-energy / (BoltzmannConstant * T)) for energy in energies) for T in temperatures]
    heat_capacity = [sum((energy - U)^2 * exp(-energy / (BoltzmannConstant * T)) for energy in energies) / (BoltzmannConstant * T^2 * sum(exp(-energy / (BoltzmannConstant * T)) for energy in energies)) for (T, U) in zip(temperatures, internal_energies)]

    # Add the lowest energy to the internal energy
    internal_energies = internal_energies .+ minimum_energy

    # Print the temperatures and internal energies in a table
    println("Temperature (K) | Internal Energy (eV) | Heat Capacity (eV/K)")
    for (T, U, C) in zip(temperatures, internal_energies, heat_capacity)
        println("$T $U $C")
    end
end

# Neighbors tests

function test_neighbors()
    # Test square lattice
    dims = (4, 4)
    println("test square lattice")
    println(sort(neighbors(:square, dims, 1, 1)) == [(1, 2), (1, 4), (2, 1), (4, 1)])  # Test neighbors of corner site
    println(sort(neighbors(:square, dims, 1, 2)) == [(1, 1), (1, 3), (2, 2), (4, 2)])  # Test neighbors of edge site
    println(sort(neighbors(:square, dims, 2, 2)) == [(1, 2), (2, 1), (2, 3), (3, 2)])  # Test neighbors of central site
    println()

    # Test hexagonal lattice
    dims = (4, 4)
    println("test hexagonal lattice")
    println("test corner sites")
    println(sort(neighbors(:hexagonal, dims, 1, 1)) == [(1, 2), (1, 4), (2, 1), (2, 2), (4, 1), (4, 2)])
    println(sort(neighbors(:hexagonal, dims, 1, 4)) == [(1, 1), (1, 3), (2, 1), (2, 4), (4, 1), (4, 4)])
    println(sort(neighbors(:hexagonal, dims, 4, 1)) == [(1, 1), (1, 4), (3, 1), (3, 4), (4, 2), (4, 4)])
    println(sort(neighbors(:hexagonal, dims, 4, 4)) == [(1, 3), (1, 4), (3, 3), (3, 4), (4, 1), (4, 3)])
    println()

    println("test edge sites")
    println(sort(neighbors(:hexagonal, dims, 1, 2)) == [(1, 1), (1, 3), (2, 2), (2, 3), (4, 2), (4, 3)])
    println(sort(neighbors(:hexagonal, dims, 2, 1)) == [(1, 1), (1, 4), (2, 2), (2, 4), (3, 1), (3, 4)])
    println(sort(neighbors(:hexagonal, dims, 1, 3)) == [(1, 2), (1, 4), (2, 3), (2, 4), (4, 3), (4, 4)])
    println(sort(neighbors(:hexagonal, dims, 3, 1)) == [(2, 1), (2, 2), (3, 2), (3, 4), (4, 1), (4, 2)])
    println(sort(neighbors(:hexagonal, dims, 2, 4)) == [(1, 3), (1, 4), (2, 1), (2, 3), (3, 3), (3, 4)])
    println(sort(neighbors(:hexagonal, dims, 4, 2)) == [(1, 1), (1, 2), (3, 1), (3, 2), (4, 1), (4, 3)])
    println(sort(neighbors(:hexagonal, dims, 3, 4)) == [(2, 1), (2, 4), (3, 1), (3, 3), (4, 1), (4, 4)])
    println(sort(neighbors(:hexagonal, dims, 4, 3)) == [(1, 2), (1, 3), (3, 2), (3, 3), (4, 2), (4, 4)])
    println()

    println("test central sites")
    println(sort(neighbors(:hexagonal, dims, 2, 2)) == [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)])
    println(sort(neighbors(:hexagonal, dims, 3, 3)) == [(2, 3), (2, 4), (3, 2), (3, 4), (4, 3), (4, 4)])
    println(sort(neighbors(:hexagonal, dims, 2, 3)) == [(1, 2), (1, 3), (2, 2), (2, 4), (3, 2), (3, 3)])
    println(sort(neighbors(:hexagonal, dims, 3, 2)) == [(2, 2), (2, 3), (3, 1), (3, 3), (4, 2), (4, 3)])
end

end # module AbstractWalkers