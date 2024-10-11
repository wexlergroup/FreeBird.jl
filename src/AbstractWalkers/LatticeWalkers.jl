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

function compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, 
                           positions::Matrix{Float64}, 
                           periodicity::Tuple{Bool, Bool, Bool}, 
                           cutoff_radii::Vector{Float64}
                           )
                           
    neighbors = Vector{Vector{Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)
    
    # Compute reciprocal lattice vectors for minimum image convention
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    a3 = supercell_lattice_vectors[:, 3]
    reciprocal_lattice_vectors = inv([a1 a2 a3])

    layers_of_neighbors = length(cutoff_radii)

    for i in 1:num_atoms
        nth_neighbors = Vector{Int}[]
        for _ in 1:layers_of_neighbors
            push!(nth_neighbors, Int[])
        end
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

                for i in 1:layers_of_neighbors
                    if distance <= cutoff_radii[i]
                        push!(nth_neighbors[i], j)
                        break
                    end 
                end


            end
        end
        
        neighbors[i] = nth_neighbors
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

    a1, a2, a3 = [lattice_vectors[:, i] for i in 1:3]

    positions = zeros(Float64, num_supercell_sites, 3)

    index = 1

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
    abstract type LatticeGeometry

The `LatticeGeometry` abstract type represents the geometry of a lattice. It has the following subtypes:

- `SquareLattice`: A square lattice.
- `TriangularLattice`: A triangular lattice.
- `GenericLattice`: A generic lattice. Currently used for non-square and non-triangular lattices.
"""
abstract type LatticeGeometry end

abstract type SquareLattice <: LatticeGeometry end

abstract type TriangularLattice <: LatticeGeometry end

abstract type GenericLattice <: LatticeGeometry end



"""
mutable struct LatticeSystem{G}

The `LatticeSystem{G}` struct represents a 3D lattice system.

# Type parameters
- `G`: The lattice geometry type. Must be a subtype of `LatticeGeometry`. See [`LatticeGeometry`](@ref).

# Fields
- `lattice_vectors::Matrix{Float64}`: The lattice vectors of the system.
- `positions::Matrix{Float64}`: The positions of the atoms in the system.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `occupations::Vector{Bool}`: A vector of booleans indicating whether each site is occupied.
- `neighbors::Vector{Tuple{Vector{Int}, Vector{Int}}}`: A vector of tuples containing the indices of the first and second nearest neighbors for each atom.
- `adsorptions::Vector{Bool}`: A vector of booleans indicating whether each site is adsorbed.

# Constructors
## Inner constructor
```julia
LatticeSystem{G}(lattice_vectors::Matrix{Float64}, 
                basis::Vector{Tuple{Float64, Float64, Float64}}, 
                supercell_dimensions::Tuple{Int64, Int64, Int64}, 
                occupations::Vector{Bool}, adsorptions::Vector{Bool}, 
                cutoff_radii::Vector{Float64},
                periodicity::Vector{Bool})
```
Create a new `LatticeSystem` with the given lattice vectors, basis, supercell dimensions, occupations, adsorptions, cutoff radii, and periodicity.

## Outer constructor for square lattice
```julia
LatticeSystem{SquareLattice}(;lattice_constant::Float64=1.0, 
                            basis=[(0.0, 0.0, 0.0)], 
                            supercell_dimensions=(4, 4, 1), 
                            periodicity=(true, true, true), 
                            cutoff_radii=[1.1, 1.5],
                            occupations=[1, 2, 3, 4],
                            adsorptions=:full)
```
Create a new `LatticeSystem` with the given square lattice parameters. The `occupations` and `adsorptions` argument 
can be a vector of integers or the symbol `:full`; the former specifies the indices of the occupied sites, while the latter specifies that all sites are occupied.

## Outer constructor for triangular lattice
```julia
LatticeSystem{TriangularLattice}(;lattice_constant::Float64=1.0, 
                                basis=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)], 
                                supercell_dimensions=(4, 2, 1), 
                                periodicity=(true, true, true), 
                                cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                occupations=[1, 2, 3, 4],
                                adsorptions=:full)
```
Create a new `LatticeSystem` with the given triangular lattice parameters. Similar to the square lattice constructor, the `occupations` and `adsorptions` argument can be a vector of integers or the symbol `:full`.

# Examples
```@repl
square_lattice  = LatticeSystem{SquareLattice}(;supercell_dimensions=(4,4,1))
triangular_lattice = LatticeSystem{TriangularLattice}(;occupations=[1,3,5,7])
```
"""
mutable struct LatticeSystem{G}
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    occupations::Vector{Bool}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}

    function LatticeSystem{G}(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        occupations::Vector{Bool},
        adsorptions::Vector{Bool},
        cutoff_radii::Vector{Float64},
    ) where G
        positions = lattice_positions(lattice_vectors, basis, supercell_dimensions)

        if length(occupations) != size(positions, 1)
            throw(ArgumentError("Length of occupations vector must match the number of lattice sites"))
        end

        if length(adsorptions) != size(positions, 1)
            throw(ArgumentError("Length of adsorptions vector must match the number of lattice sites"))
        end

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii)
        
        return new{G}(lattice_vectors, positions, basis, supercell_dimensions, periodicity, occupations, neighbors, adsorptions)
    end
end

function LatticeSystem{SquareLattice}(;lattice_constant::Float64=1.0,
                                      basis::Vector{Tuple{Float64, Float64, Float64}}=[(0.0, 0.0, 0.0)],
                                      supercell_dimensions::Tuple{Int64, Int64, Int64}=(4, 4, 1),
                                      periodicity::Tuple{Bool, Bool, Bool}=(true, true, true),
                                      cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                      occupations::Union{Vector{Int}, Symbol}=[1, 2, 3, 4],
                                      adsorptions::Union{Vector{Int}, Symbol}=:full,
                                      )

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 lattice_constant 0.0; 0.0 0.0 1.0]
    dim = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3]
    lattice_occupations = zeros(Bool, dim*length(basis))
    lattice_adsorptions = zeros(Bool, dim*length(basis))

    if occupations == :full
        lattice_occupations = [true for i in 1:dim*length(basis)]
    else
        for i in occupations
            lattice_occupations[i] = true
        end
    end

    if adsorptions == :full
        lattice_adsorptions = [true for i in 1:dim*length(basis)]
    else
        for i in adsorptions
            lattice_adsorptions[i] = true
        end
    end

    return LatticeSystem{SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_occupations, lattice_adsorptions, cutoff_radii)
end


function LatticeSystem{TriangularLattice}(;lattice_constant::Float64=1.0,
                                          basis::Vector{Tuple{Float64, Float64, Float64}}=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)],
                                          supercell_dimensions::Tuple{Int64, Int64, Int64}=(4, 2, 1),
                                          periodicity::Tuple{Bool, Bool, Bool}=(true, true, true),
                                          cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                          occupations::Union{Vector{Int}, Symbol}=[1, 2, 3, 4],
                                          adsorptions::Union{Vector{Int}, Symbol}=:full,
                                          )

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 sqrt(3)*lattice_constant 0.0; 0.0 0.0 1.0]
    dim = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3]
    lattice_occupations = zeros(Bool, dim * length(basis))
    lattice_adsorptions = zeros(Bool, dim * length(basis))

    if occupations == :full
        lattice_occupations = [true for i in 1:dim*length(basis)]
    else
        for i in occupations
            lattice_occupations[i] = true
        end
    end

    if adsorptions == :full
        lattice_adsorptions = [true for i in 1:dim*length(basis)]
    else
        for i in adsorptions
            lattice_adsorptions[i] = true
        end
    end

    return LatticeSystem{TriangularLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_occupations, lattice_adsorptions, cutoff_radii)
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

function Base.show(io::IO, walker::LatticeWalker)
    println(io, "LatticeWalker(")
    println(io, "    configuration: ", walker.configuration)
    println(io, "    energy: ", walker.energy)
    println(io, "    iter: ", walker.iter, ")")
end

function Base.show(io::IO, walker::Vector{LatticeWalker})
    println(io, "Vector{LatticeWalker}(", length(walker), "):")
    for (ind, w) in enumerate(walker)
        println(io, "[", ind, "] ", w)
    end
end
