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
lattice_positions(lattice_vectors::Matrix{Float64}, basis::Vector{Tuple{Float64, Float64, Float64}}, supercell_dimensions::Tuple{Int64, Int64, Int64})

Compute the positions of atoms in a 3D lattice.

# Arguments
- `lattice_vectors::Matrix{Float64}`: The lattice vectors of the system.
- `basis::Vector{Tuple{Float64, Float64, Float64}}`: The basis of the system.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.

# Returns
- `positions::Matrix{Float64}`: The positions of the atoms in the supercell.

"""
function lattice_positions(lattice_vectors::Matrix{Float64}, 
                           basis::Vector{Tuple{Float64, Float64, Float64}}, 
                           supercell_dimensions::Tuple{Int64, Int64, Int64},
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
    
    return positions
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
        positions = lattice_positions(lattice_vectors, basis, supercell_dimensions)

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


function Base.show(io::IO, lattice::LatticeSystem)
    println(io, "LatticeSystem:")
    println(io, "    lattice_vectors      : ", lattice.lattice_vectors)
    println(io, "    positions            : ", lattice.positions)
    println(io, "    supercell_dimensions : ", lattice.supercell_dimensions)
    println(io, "    occupations          : ", lattice.occupations)
    println(io, "    adsorptions          : ", lattice.adsorptions)
    println(io, "    neighbors            : ")
    if length(lattice.neighbors) > 10
        for i in 1:5
            println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
        end
        println(io, "        ⋮")
        for i in length(lattice.neighbors)-4:length(lattice.neighbors)
            println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
        end
    else
        for i in 1:length(lattice.neighbors)
            println(io, "        site ", i, ": ", "nearest = ", lattice.neighbors[i][1], ", next-nearest = ", lattice.neighbors[i][2])
        end
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

mutable struct LatticeWalker <: AbstractWalker
    configuration::LatticeSystem
    energy::typeof(0.0u"eV")
    iter::Int64
    function LatticeWalker(configuration::LatticeSystem; energy=0.0u"eV", iter=0)
        return new(configuration, energy, iter)
    end
end
