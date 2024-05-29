abstract type LatticeWalkers end

"""
    compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, cutoff_radii::Tuple{Float64, Float64})

Compute the nearest and next-nearest neighbors for each atom in a 2D lattice.

# Arguments
- `supercell_lattice_vectors::Matrix{Float64}`: The lattice vectors of the supercell.
- `positions::Matrix{Float64}`: The positions of the atoms in the supercell.
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.

# Returns
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.

"""

function compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, cutoff_radii::Tuple{Float64, Float64})
    neighbors = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)
    
    # Compute reciprocal lattice vectors for minimum image convention
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    reciprocal_lattice_vectors = inv([a1 a2])

    first_nearest_distance = cutoff_radii[1]
    second_nearest_distance = cutoff_radii[2]

    for i in 1:num_atoms
        first_neighbors = Int[]
        second_neighbors = Int[]
        pos_i = positions[i, :]
        
        for j in 1:num_atoms
            if i != j
                pos_j = positions[j, :]
                dx = pos_j[1] - pos_i[1]
                dy = pos_j[2] - pos_i[2]

                # Apply minimum image convention using reciprocal lattice vectors
                dr = [dx, dy]
                fractional_dr = reciprocal_lattice_vectors * dr
                fractional_dr .= fractional_dr .- round.(fractional_dr)
                dr = supercell_lattice_vectors * fractional_dr

                distance = norm(dr)
                
                if distance <= first_nearest_distance
                    push!(first_neighbors, j)
                elseif distance <= second_nearest_distance
                    push!(second_neighbors, j)
                end
            end
        end
        
        neighbors[i] = (first_neighbors, second_neighbors)
    end
    
    return neighbors
end

"""
    mutable struct LatticeSystem

The `LatticeSystem` struct represents a 3D lattice system.

# Fields
- `lattice_vectors::Matrix{Float64}`: The lattice vectors of the system.
- `positions::Matrix{Float64}`: The positions of the atoms in the system.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `occupations::Vector{Bool}`: A vector of booleans indicating whether each site is occupied.
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.

# Constructors
```julia
LatticeSystem(lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64}, occupations::Vector{Bool}, cutoff_radii::Tuple{Float64, Float64})
```
Create a new `LatticeSystem` with the given lattice vectors, basis, supercell dimensions, occupations, and cutoff radii.

"""

mutable struct LatticeSystem
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    supercell_dimensions::Tuple{Int64, Int64}
    occupations::Vector{Bool}
    neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}

    function LatticeSystem(lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64}, occupations::Vector{Bool}, cutoff_radii::Tuple{Float64, Float64})
        num_basis_sites = length(basis)
        num_supercell_sites = supercell_dimensions[1] * supercell_dimensions[2] * num_basis_sites
        
        positions = zeros(Float64, num_supercell_sites, 2)
        index = 1

        a1 = lattice_vectors[:, 1]
        a2 = lattice_vectors[:, 2]

        for j in 1:supercell_dimensions[2]
            for i in 1:supercell_dimensions[1]
                for (bx, by) in basis
                    x = (i - 1) * a1[1] + (j - 1) * a2[1] + bx
                    y = (i - 1) * a1[2] + (j - 1) * a2[2] + by
                    positions[index, :] = [x, y]
                    index += 1
                end
            end
        end

        if length(occupations) != size(positions, 1)
            throw(ArgumentError("Length of occupations vector must match the number of lattice sites"))
        end

        neighbors = compute_neighbors(lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2]]), positions, cutoff_radii)
        
        return new(lattice_vectors, positions, supercell_dimensions, occupations, neighbors)
    end
end

"""
    mutable struct LatticeWalker

The `LatticeWalker` struct represents a walker on a 3D lattice.

# Fields
- `configuration::LatticeSystem`: The configuration of the walker.
- `energy::Float64`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.

# Constructor
```julia
LatticeWalker(configuration::LatticeSystem; energy=0.0, iter=0)
```
Create a new `LatticeWalker` with the given configuration and optional energy and iteration number.

"""

mutable struct LatticeWalker
    configuration::LatticeSystem
    energy::Float64
    iter::Int64
    function LatticeWalker(configuration::LatticeSystem; energy=0.0, iter=0)
        return new(configuration, energy, iter)
    end
end

function nearest_neighbors(lattice_type::Symbol, dims::Tuple{Int64, Int64}, i::Int64, j::Int64)
    if lattice_type == :square
        return [(mod(i + dx - 1, dims[1]) + 1, mod(j + dy - 1, dims[2]) + 1) 
                for (dx, dy) in [(0, 1), (0, -1), (1, 0), (-1, 0)]]
    elseif lattice_type == :hexagonal
        even_row_neighbors = [(0, -1), (0, 1), (-1, -1), (-1, 0), (1, -1), (1, 0)]
        odd_row_neighbors = [(0, -1), (0, 1), (-1, 0), (-1, 1), (1, 0), (1, 1)]
        offsets = iseven(i) ? even_row_neighbors : odd_row_neighbors
        return [(mod(i + dx - 1, dims[1]) + 1, mod(j + dy - 1, dims[2]) + 1) for (dx, dy) in offsets]
    else
        error("Lattice type not implemented")
    end
end

function next_nearest_neighbors(lattice_type::Symbol, dims::Tuple{Int64, Int64}, i::Int64, j::Int64)
    if lattice_type == :square
        return [(mod(i + dx - 1, dims[1]) + 1, mod(j + dy - 1, dims[2]) + 1)
                for (dx, dy) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]]
    elseif lattice_type == :hexagonal
        # TODO: Implement hexagonal lattice next-nearest-neighbor interaction
        error("Hexagonal lattice next-nearest-neighbor interaction not implemented yet")
    else
        error("Lattice type not implemented")
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
                for (i2, j2) in nearest_neighbors(at.lattice_type, at.dimensions, i, j)
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
                for (i2, j2) in next_nearest_neighbors(at.lattice_type, at.dimensions, i, j)
                    e_nnn_interaction += at.site_occupancy[i2, j2] ? lg.nnn_interaction_energy / 2 : 0.0u"eV"
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

    # Generate a list of energies and configurations
    energies = [walker.energy for walker in lg_walkers.walkers]
    configurations = [walker.configuration for walker in lg_walkers.walkers]

    return energies, configurations
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