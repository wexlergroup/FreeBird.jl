abstract type LatticeWalkers end

"""
compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, cutoff_radii::Tuple{Float64, Float64}, periodicity::Vector{Bool})

Compute the nearest and next-nearest neighbors for each atom in a 3D lattice.

# Arguments
- `supercell_lattice_vectors::Matrix{Float64}`: The lattice vectors of the supercell.
- `positions::Matrix{Float64}`: The positions of the atoms in the supercell.
- `periodicity::Vector{Bool}`: A Boolean vector of length three indicating periodicity in each dimension (true for periodic, false for non-periodic).
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.

# Returns
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.

"""

function compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, positions::Matrix{Float64}, periodicity::Tuple{Bool, Bool, Bool}, cutoff_radii::Tuple{Float64, Float64})
    neighbors = Vector{Tuple{Vector{Int}, Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)
    
    # Compute reciprocal lattice vectors for minimum image convention
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    a3 = supercell_lattice_vectors[:, 3]
    reciprocal_lattice_vectors = inv([a1 a2 a3])

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
                dz = pos_j[3] - pos_i[3]

                # Apply minimum image convention using reciprocal lattice vectors
                dr = [dx, dy, dz]
                fractional_dr = reciprocal_lattice_vectors * dr
                
                for k in 1:3
                    if periodicity[k]
                        fractional_dr[k] -= round(fractional_dr[k])
                    end
                end
                
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
- `adsorptions::Vector{Bool}`: A vector of booleans indicating whether each site is adsorbed.

# Constructors
```julia
LatticeSystem(lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64, Int64}, occupations::Vector{Bool}, adsorptions::Vector{Bool}, cutoff_radii::Tuple{Float64, Float64}, periodicity::Vector{Bool})
```
Create a new `LatticeSystem` with the given lattice vectors, basis, supercell dimensions, occupations, adsorptions, cutoff radii, and periodicity.

"""

mutable struct LatticeSystem
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    occupations::Vector{Bool}
    neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}
    adsorptions::Vector{Bool}

    function LatticeSystem(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        occupations::Vector{Bool},
        adsorptions::Vector{Bool},
        cutoff_radii::Tuple{Float64, Float64}
    )
        num_basis_sites = length(basis)
        num_supercell_sites = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3] * num_basis_sites
        
        positions = zeros(Float64, num_supercell_sites, 3)
        index = 1

        a1 = lattice_vectors[:, 1]
        a2 = lattice_vectors[:, 2]
        a3 = lattice_vectors[:, 3]

        for k in 1:supercell_dimensions[3]
            for j in 1:supercell_dimensions[2]
                for i in 1:supercell_dimensions[1]
                    for (bx, by, bz) in basis
                        x = (i - 1) * a1[1] + (j - 1) * a2[1] + (k - 1) * a3[1] + bx
                        y = (i - 1) * a1[2] + (j - 1) * a2[2] + (k - 1) * a3[2] + by
                        z = (i - 1) * a1[3] + (j - 1) * a2[3] + (k - 1) * a3[3] + bz
                        positions[index, :] = [x, y, z]
                        index += 1
                    end
                end
            end
        end

        if length(occupations) != size(positions, 1)
            throw(ArgumentError("Length of occupations vector must match the number of lattice sites"))
        end

        if length(adsorptions) != size(positions, 1)
            throw(ArgumentError("Length of adsorptions vector must match the number of lattice sites"))
        end

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii)
        
        return new(lattice_vectors, positions, supercell_dimensions, occupations, neighbors, adsorptions)
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

"""
    interaction_energy(at::LatticeSystem, lg::LGHamiltonian)

Compute the interaction energy of a lattice configuration using the Hamiltonian parameters.

# Arguments
- `at::LatticeSystem`: The lattice configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.

# Returns
- `e_interaction::Float64`: The interaction energy of the lattice configuration.

"""

function interaction_energy(at::LatticeSystem, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64)
    e_adsorption = sum(at.occupations .& at.adsorptions) * adsorption_energy
    e_nn = 0.0
    e_nnn = 0.0

    for index in 1:length(at.occupations)
        if at.occupations[index]
            # Compute nearest-neighbor interaction energy
            for nn in at.neighbors[index][1]
                if at.occupations[nn]
                    e_nn += nn_energy / 2
                end
            end

            # Compute next-nearest-neighbor interaction energy
            for nnn in at.neighbors[index][2]
                if at.occupations[nnn]
                    e_nnn += nnn_energy / 2
                end
            end
        end
    end

    e_interaction = e_adsorption + e_nn + e_nnn
    return e_interaction
end

# """
#     assign_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)

# Assigns energies to a list of `Lattice2DWalker` objects based on the Hamiltonian parameters.

# # Arguments
# - `walkers::Vector{Lattice2DWalker}`: A list of `Lattice2DWalker` objects to assign energies to.
# - `lg::LGHamiltonian`: The Hamiltonian parameters used for energy calculation.

# # Returns
# - `walkers::Vector{Lattice2DWalker}`: The updated list of `Lattice2DWalker` objects with assigned energies.

# """

# function assign_energies!(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
#     for walker in walkers
#         e_adsorption = walker.configuration.num_occ_sites * lg.adsorption_energy
#         e_interaction = interaction_energy(walker.configuration, lg)  # nearest-neighbor and next-nearest-neighbor interactions
#         e_total = e_adsorption + e_interaction
#         walker.energy = e_total
#     end
#     return walkers
# end

# """
#     struct Lattice2DWalkers <: LatticeWalkers

# The `Lattice2DWalkers` struct contains a list of `Lattice2DWalker` objects and the Hamiltonian parameters
#     that defines the interactions between the particles in the walkers.

# # Fields
# - `walkers::Vector{Lattice2DWalker}`: The list of `Lattice2DWalker` objects.
# - `lg::LGHamiltonian`: The Hamiltonian parameters.

# # Constructors
# - `Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)`: Constructs a new `Lattice2DWalkers`
#     object with the given walkers and Hamiltonian parameters. The energies of the walkers are automatically
#     assigned using the Hamiltonian parameters.

# """

# struct Lattice2DWalkers <: LatticeWalkers
#     walkers::Vector{Lattice2DWalker}
#     lg::LGHamiltonian
#     function Lattice2DWalkers(walkers::Vector{Lattice2DWalker}, lg::LGHamiltonian)
#         assign_energies!(walkers, lg)
#         return new(walkers, lg)
#     end
# end

"""
    exact_enumeration(primitive_lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64}, number_occupied_sites::Int64, adsorption_energy::Float64, nn_energy::Float64, nnn_energy::Float64, cutoff_radii::Tuple{Float64, Float64})

Enumerate all possible configurations of a lattice system and compute the energy of each configuration.

# Arguments
- `primitive_lattice_vectors::Matrix{Float64}`: The primitive lattice vectors of the system.
- `basis::Vector{Tuple{Float64, Float64}}`: The basis of the system.
- `supercell_dimensions::Tuple{Int64, Int64}`: The dimensions of the supercell.
- `number_occupied_sites::Int64`: The number of occupied sites in each configuration.
- `adsorption_energy::Float64`: The adsorption energy of the particles.
- `nn_energy::Float64`: The nearest-neighbor interaction energy.
- `nnn_energy::Float64`: The next-nearest-neighbor interaction energy.
- `cutoff_radii::Tuple{Float64, Float64}`: The cutoff radii for the first and second nearest neighbors.

# Returns
- `energies::Vector{Float64}`: A vector of the energies of each configuration.
- `configurations::Vector{Vector{Bool}}`: A vector of the configurations of the lattice system.

"""

function exact_enumeration(
    primitive_lattice_vectors::Matrix{Float64},
    basis::Vector{Tuple{Float64, Float64, Float64}},
    supercell_dimensions::Tuple{Int64, Int64, Int64},
    periodicity::Tuple{Bool, Bool, Bool},
    number_occupied_sites::Int64,
    adsorptions::Vector{Bool},
    adsorption_energy::Float64,
    nn_energy::Float64,
    nnn_energy::Float64,
    cutoff_radii::Tuple{Float64, Float64}
)
    K, L, M = supercell_dimensions
    num_basis_sites = length(basis)
    total_sites = K * L * M * num_basis_sites

    # Generate all possible occupation configurations
    all_configs = collect(combinations(1:total_sites, number_occupied_sites))

    # Generate occupation vectors from configurations
    all_occupation_vectors = Vector{Vector{Bool}}()
    for config in all_configs
        occupations = falses(total_sites)
        occupations[config] .= true
        push!(all_occupation_vectors, Vector{Bool}(occupations))  # Convert BitVector to Vector{Bool}
    end

    # Generate LatticeSystem objects for each configuration
    lattices = [LatticeSystem(primitive_lattice_vectors, basis, supercell_dimensions, periodicity, occupations, adsorptions, cutoff_radii) for occupations in all_occupation_vectors]

    # Generate LatticeWalker objects for each lattice system
    walkers = [LatticeWalker(lattice) for lattice in lattices]

    # Compute energies for each walker
    for walker in walkers
        e_interaction = interaction_energy(walker.configuration, adsorption_energy, nn_energy, nnn_energy)
        walker.energy = e_interaction
    end

    # Extract energies and configurations
    energies = [walker.energy for walker in walkers]
    configurations = [walker.configuration for walker in walkers]

    return energies, configurations
end

# function compute_internal_energy_versus_temperature(L::Int64, N::Int64, T_min::typeof(1.0u"K"), T_max::typeof(100.0u"K"), num_points::Int64, lattice_type::Symbol=:square)
#     # Generate all possible configurations with 8 occupied sites
#     energies = exact_enumeration(L, L, N, lattice_type, LGHamiltonian(-0.04136319965u"eV", -0.01034079991u"eV", -(15 / 64) * 0.01034079991u"eV"))

#     # Compute energy relative to the lowest energy
#     minimum_energy = minimum(energies)
#     energies = energies .- minimum_energy

#     # Compute the partition function
#     BoltzmannConstant = 8.617_333_262e-5u"eV/K"
#     temperatures = range(T_min, T_max, length=num_points)
#     internal_energies = [sum(energy * exp(-energy / (BoltzmannConstant * T)) for energy in energies) / sum(exp(-energy / (BoltzmannConstant * T)) for energy in energies) for T in temperatures]
#     heat_capacity = [sum((energy - U)^2 * exp(-energy / (BoltzmannConstant * T)) for energy in energies) / (BoltzmannConstant * T^2 * sum(exp(-energy / (BoltzmannConstant * T)) for energy in energies)) for (T, U) in zip(temperatures, internal_energies)]

#     # Add the lowest energy to the internal energy
#     internal_energies = internal_energies .+ minimum_energy

#     # Print the temperatures and internal energies in a table
#     println("Temperature (K) | Internal Energy (eV) | Heat Capacity (eV/K)")
#     for (T, U, C) in zip(temperatures, internal_energies, heat_capacity)
#         println("$T $U $C")
#     end
# end