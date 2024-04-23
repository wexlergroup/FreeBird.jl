abstract type LatticeWalkers end

"""
    mutable struct Lattice2DSystem

The `Lattice2DSystem` struct represents a 2D lattice system.

# Fields
- `lattice_type::Symbol`: The type of lattice (e.g., :square, :hexagonal).
- `dimensions::Tuple{Int64, Int64}`: The dimensions of the lattice (rows x columns for 2D).
- `num_occ_sites::Int64`: The number of occupied sites.
- `site_occupancy::Matrix{Bool}`: A matrix to track whether a lattice site is occupied.

# Constructors
```julia
Lattice2DSystem(lattice_type::Symbol, site_occupancy::Matrix{Bool})
```
Create a new `Lattice2DSystem` with the given lattice type and site occupancy matrix.

```julia
Lattice2DSystem(lattice_type::Symbol, dims::Tuple{Int64, Int64}, num_occ_sites::Int64)
```
Create a new `Lattice2DSystem` with the given lattice type, dimensions, and number of occupied sites.

"""

mutable struct Lattice2DSystem
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

"""
    mutable struct Lattice2DWalker

The `Lattice2DWalker` struct represents a walker on a 2D lattice.

# Fields
- `configuration::Lattice2DSystem`: The configuration of the walker.
- `energy::typeof(0.0u"eV")`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.

# Constructor
```julia
Lattice2DWalker(configuration::Lattice2DSystem; energy=0.0u"eV", iter=0)
```
Create a new `Lattice2DWalker` with the given configuration and optional energy and iteration number.

"""

mutable struct Lattice2DWalker
    configuration::Lattice2DSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    function Lattice2DWalker(configuration::Lattice2DSystem; energy=0.0u"eV", iter=0)
        return new(configuration, energy, iter)
    end
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

"""
    assign_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)

Assigns energies to a list of `Lattice2DWalker` objects based on the Hamiltonian parameters.

# Arguments
- `walkers::Vector{Lattice2DWalker}`: A list of `Lattice2DWalker` objects to assign energies to.
- `lg::LGHamiltonian`: The Hamiltonian parameters used for energy calculation.

# Returns
- `walkers::Vector{Lattice2DWalker}`: The updated list of `Lattice2DWalker` objects with assigned energies.

"""

function assign_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
    for walker in walkers
        e_adsorption = walker.configuration.num_occ_sites * lg.adsorption_energy
        e_interaction = interaction_energy(walker.configuration, lg)  # nearest-neighbor and next-nearest-neighbor interactions
        e_total = e_adsorption + e_interaction
        walker.energy = e_total
    end
    return walkers
end

"""
    struct Lattice2DWalkers <: LatticeWalkers

The `Lattice2DWalkers` struct contains a list of `Lattice2DWalker` objects and the Hamiltonian parameters
    that defines the interactions between the particles in the walkers.

# Fields
- `walkers::Vector{Lattice2DWalker}`: The list of `Lattice2DWalker` objects.
- `lg::LGHamiltonian`: The Hamiltonian parameters.

# Constructors
- `Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)`: Constructs a new `Lattice2DWalkers`
    object with the given walkers and Hamiltonian parameters. The energies of the walkers are automatically
    assigned using the Hamiltonian parameters.

"""

struct Lattice2DWalkers <: LatticeWalkers
    walkers::Vector{Lattice2DWalker}
    lg::LGHamiltonian
    function Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
        assign_energies!(walkers, lg)
        return new(walkers, lg)
    end
end

function exact_enumeration(L::Int64, M::Int64, N::Int64, lattice_type::Symbol, lg::LGHamiltonian)
    # Generate a list of all sites
    sites = [(i, j) for i in 1:L for j in 1:M]

    # Generate all combinations of N sites
    all_configs = collect(combinations(sites, N))

    # Convert each configuration to a site occupancy matrix
    all_configs = [reshape([in((i, j), config) for i in 1:L, j in 1:M], L, M) for config in all_configs]

    # Generate Lattice2DSystem objects for each configuration
    lattices = [Lattice2DSystem(lattice_type, config) for config in all_configs]

    # Compute the energy of each configuration
    if lattice_type == :hexagonal && lg.nnn_interaction_energy != 0.0u"eV"
        error("Next-nearest-neighbor interaction energy not implemented for hexagonal lattice")
    end
    walkers = [Lattice2DWalker(lattice) for lattice in lattices]
    lg_walkers = Lattice2DWalkers(walkers, lg)

    # Generate a list of energies
    energies = [walker.energy for walker in lg_walkers.walkers]

    return energies
end

function compute_internal_energy_versus_temperature(L::Int64, N::Int64, T_min::typeof(1.0u"K"), T_max::typeof(100.0u"K"), num_points::Int64, lattice_type::Symbol=:square)
    # Generate all possible configurations with 8 occupied sites
    energies = exact_enumeration(L, L, N, lattice_type, LGHamiltonian(-0.04136319965u"eV", -0.01034079991u"eV", -(15 / 64) * 0.01034079991u"eV"))

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