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

abstract type AbstractLattice end

"""
mutable struct SLattice{G}

The `SLattice{G}` struct represents a 3D lattice system.

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
SLattice{G}(lattice_vectors::Matrix{Float64}, 
                basis::Vector{Tuple{Float64, Float64, Float64}}, 
                supercell_dimensions::Tuple{Int64, Int64, Int64}, 
                occupations::Vector{Bool}, adsorptions::Vector{Bool}, 
                cutoff_radii::Vector{Float64},
                periodicity::Vector{Bool})
```
Create a new `SLattice` with the given lattice vectors, basis, supercell dimensions, occupations, adsorptions, cutoff radii, and periodicity.

## Outer constructor for square lattice
```julia
SLattice{SquareLattice}(;lattice_constant::Float64=1.0, 
                            basis=[(0.0, 0.0, 0.0)], 
                            supercell_dimensions=(4, 4, 1), 
                            periodicity=(true, true, true), 
                            cutoff_radii=[1.1, 1.5],
                            occupations=[1, 2, 3, 4],
                            adsorptions=:full)
```
Create a new `SLattice` with the given square lattice parameters. The `occupations` and `adsorptions` argument 
can be a vector of integers or the symbol `:full`; the former specifies the indices of the occupied sites, while the latter specifies that all sites are occupied.

## Outer constructor for triangular lattice
```julia
SLattice{TriangularLattice}(;lattice_constant::Float64=1.0, 
                                basis=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)], 
                                supercell_dimensions=(4, 2, 1), 
                                periodicity=(true, true, true), 
                                cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                occupations=[1, 2, 3, 4],
                                adsorptions=:full)
```
Create a new `SLattice` with the given triangular lattice parameters. Similar to the square lattice constructor, the `occupations` and `adsorptions` argument can be a vector of integers or the symbol `:full`.

# Examples
```@repl
square_lattice  = SLattice{SquareLattice}(;supercell_dimensions=(4,4,1))
triangular_lattice = SLattice{TriangularLattice}(;occupations=[1,3,5,7])
```
"""
mutable struct SLattice{G} <: AbstractLattice
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    occupations::Vector{Bool}
    cutoff_radii::Vector{Float64}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}

    function SLattice{G}(
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
        
        return new{G}(lattice_vectors, positions, basis, supercell_dimensions, periodicity, occupations, cutoff_radii, neighbors, adsorptions)
    end
end

function SLattice{SquareLattice}(;lattice_constant::Float64=1.0,
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

    return SLattice{SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_occupations, lattice_adsorptions, cutoff_radii)
end


function SLattice{TriangularLattice}(;lattice_constant::Float64=1.0,
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

    return SLattice{TriangularLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_occupations, lattice_adsorptions, cutoff_radii)
end








"""
    mutable struct MLattice{C,G}

A mutable struct representing a lattice with the following fields:

- `lattice_vectors::Matrix{Float64}`: The lattice vectors defining the unit cell.
- `positions::Matrix{Float64}`: The positions of the lattice points.
- `basis::Vector{Tuple{Float64, Float64, Float64}}`: The basis vectors within the unit cell.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `periodicity::Tuple{Bool, Bool, Bool}`: The periodicity in each dimension.
- `components::Vector{Vector{Bool}}`: The components of the lattice.
- `neighbors::Vector{Vector{Vector{Int}}}`: The neighbors of each lattice point.
- `adsorptions::Vector{Bool}`: The adsorption sites on the lattice.

# Inner Constructor

    MLattice{C,G}(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        components::Vector{Vector{Bool}},
        adsorptions::Vector{Bool},
        cutoff_radii::Vector{Float64},
    ) where {C,G}

Creates an `MLattice` instance with the specified parameters. The constructor performs the following steps:

1. Validates that the number of components matches the expected value `C`.
2. Computes the positions of the lattice points using `lattice_positions`.
3. Computes the supercell lattice vectors.
4. Computes the neighbors of each lattice point using `compute_neighbors`.

Throws an `ArgumentError` if the number of components does not match `C`.

# Outer Constructors

    MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
                               basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
                               supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
                               periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                               cutoff_radii::Vector{Float64}=[1.1, 1.5],
                               components::Union{Vector{Int},Symbol}=:equal,
                               adsorptions::Union{Vector{Int},Symbol}=:full)

Constructs a square lattice with the specified parameters. The `components` and `adsorptions` arguments can be a vector of integers specifying
the indices of the occupied sites, or a symbol. If `components` is `:equal`, the lattice is divided into `C` equal components when possible, or 
nearest to equal components otherwise. If `adsorptions` is `:full`, all sites are classified as adsorption sites.

## Returns
- `MLattice{C,SquareLattice}`: A square lattice object with `C` components.

"""
mutable struct MLattice{C,G} <: AbstractLattice
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    components::Vector{Vector{Bool}}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}

    function MLattice{C,G}(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        components::Vector{Vector{Bool}},
        adsorptions::Vector{Bool},
        cutoff_radii::Vector{Float64},
    ) where {C,G}

        num_components = length(components)

        if num_components != C
            throw(ArgumentError("For a $C-component system, got $num_components components!"))
        end

        positions = lattice_positions(lattice_vectors, basis, supercell_dimensions)

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii)
        
        return new{C,G}(lattice_vectors, positions, basis, supercell_dimensions, periodicity, components, neighbors, adsorptions)
    end
end

function split_into_subarrays(arr::AbstractVector, N::Int)
    n = length(arr)  # Total number of elements
    base_size = div(n, N)  # Base size of each subarray
    remainder = mod(n, N)  # Remaining elements to distribute

    subarrays = Vector{Vector{eltype(arr)}}()
    idx = 1

    for i in 1:N
        # Determine the size of the current subarray
        current_size = base_size + (i <= remainder ? 1 : 0)
        push!(subarrays, arr[idx:idx + current_size - 1])
        idx += current_size
    end

    return subarrays
end

function MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
                                    basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
                                    supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
                                    periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                    cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                    components::Union{Vector{Int},Symbol}=:equal,
                                    adsorptions::Union{Vector{Int},Symbol}=:full,
                                ) where C

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 lattice_constant 0.0; 0.0 0.0 1.0]
    dim = prod(supercell_dimensions)
    lattice_adsorptions = zeros(Bool, dim * length(basis))

    if adsorptions == :full
        lattice_adsorptions = [true for i in 1:dim*length(basis)]
    else
        for i in adsorptions
            lattice_adsorptions[i] = true
        end
    end

    
    if components == :equal
        lattice_comp = Vector{Vector{Bool}}(undef, C)
        comps = split_into_subarrays(1:dim, C)
        for i in 1:C
            lattice_comp[i] = [false for i in 1:dim]
            for j in comps[i]
                lattice_comp[i][j] = true
            end
        end
    else
        lattice_comp = components
    end



    return MLattice{C,SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, lattice_comp, lattice_adsorptions, cutoff_radii)
end

"""
    mutable struct LatticeWalker

The `LatticeWalker` struct represents a walker on a 3D lattice.

# Fields
- `configuration::AbstractLattice`: The configuration of the walker.
- `energy::Float64`: The energy of the walker.
- `iter::Int64`: The current iteration number of the walker.

# Constructor
```julia
LatticeWalker(configuration::AbstractLattice; energy=0.0, iter=0)
```
Create a new `LatticeWalker` with the given configuration and optional energy and iteration number.

"""


number_of_lattice_components(lattice::SLattice) = 1

number_of_lattice_components(lattice::MLattice{C,G}) where {C,G} = C

mutable struct LatticeWalker{C} <: AbstractWalker
    configuration::AbstractLattice
    energy::typeof(0.0u"eV")
    iter::Int64
    function LatticeWalker(configuration::AbstractLattice; energy=0.0u"eV", iter=0)
        num_comp = number_of_lattice_components(configuration)
        return new{num_comp}(configuration, energy, iter)
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