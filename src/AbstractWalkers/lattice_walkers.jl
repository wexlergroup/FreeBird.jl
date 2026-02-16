"""
compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, 
                  positions::Matrix{Float64}, 
                  cutoff_radii::Vector{Float64}, 
                  periodicity::Tuple{Bool, Bool, Bool})

Compute the nearest and next-nearest neighbors for each atom in a 3D lattice.

# Arguments
- `supercell_lattice_vectors::Matrix{Float64}`: The lattice vectors of the supercell.
- `positions::Matrix{Float64}`: The positions of the atoms in the supercell.
- `periodicity::Tuple{Bool, Bool, Bool}`: A Boolean tuple of length three indicating periodicity in each dimension (true for periodic, false for non-periodic).
- `cutoff_radii::Vector{Float64}`: The cutoff radii for the *index*-th nearest neighbors.

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

function get_lattice_positions(lattice_vectors::Matrix{Float64}, supercell_dimensions::Tuple{Int64, Int64, Int64})
    num_supercell_sites = supercell_dimensions[1] * supercell_dimensions[2] * supercell_dimensions[3]

    a1, a2, a3 = [lattice_vectors[:, i] for i in 1:3]

    positions = zeros(Float64, num_supercell_sites, 3)

    index = 1

    for k in 1:supercell_dimensions[3]
        for j in 1:supercell_dimensions[2]
            for i in 1:supercell_dimensions[1]
                x = (i - 1) * a1[1] + (j - 1) * a2[1] + (k - 1) * a3[1]
                y = (i - 1) * a1[2] + (j - 1) * a2[2] + (k - 1) * a3[2]
                z = (i - 1) * a1[3] + (j - 1) * a2[3] + (k - 1) * a3[3]
                positions[index, :] = [x, y, z]
                index += 1
            end
        end
    end
    
    return positions
end

function find_n_cutoff_radii(positions::Matrix{Float64}, num_nearest_neighbors::Int64)
    x_dist = positions[2, 1] - positions[1, 1]
    y_dist = positions[2, 2] - positions[1, 2]
    first_nn = sqrt(x_dist^2 + y_dist^2)
    cutoff_radii = vcat([first_nn], zeros(Float64, num_nearest_neighbors - 1))
    for i in 2:num_nearest_neighbors
        if iseven(i)
            radii = (i / 2) * sqrt(2) * first_nn
        else
            radii = ((i + 1) / 2) * first_nn
        end
        cutoff_radii[i] = radii
    end
    return cutoff_radii
end

function compute_neighbors(supercell_lattice_vectors::Matrix{Float64}, 
                           positions::Matrix{Float64}, 
                           periodicity::Tuple{Bool, Bool, Bool}, 
                           cutoff_radii::Vector{Float64})
    neighbors = Vector{Vector{Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)
    
    # Extract lattice vectors
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    a3 = supercell_lattice_vectors[:, 3]
    
    # Handle 2D case: only invert the non-zero part
    if !periodicity[3] && all(a3 .== 0)
        # 2D system - construct 2x2 inverse for x,y only
        lattice_2d = [a1[1:2] a2[1:2]]
        inv_lattice_2d = inv(lattice_2d)
        reciprocal_lattice_vectors = zeros(3, 3)
        reciprocal_lattice_vectors[1:2, 1:2] = inv_lattice_2d
    else
        # 3D system
        reciprocal_lattice_vectors = inv([a1 a2 a3])
    end
    
    layers_of_neighbors = length(cutoff_radii)
    
    for i in 1:num_atoms
        nth_neighbors = [Int[] for _ in 1:layers_of_neighbors]
        pos_i = positions[i, :]
        
        for j in 1:num_atoms
            if i != j
                pos_j = positions[j, :]
                dr = pos_j - pos_i
                
                # Apply minimum image convention
                if periodicity[1] || periodicity[2]
                    fractional_dr = reciprocal_lattice_vectors * dr
                    for k in 1:3
                        if periodicity[k]
                            fractional_dr[k] -= round(fractional_dr[k])
                        end
                    end
                    dr = supercell_lattice_vectors * fractional_dr
                end
                
                distance = norm(dr)
                
                # Assign to neighbor shell
                for layer in 1:layers_of_neighbors
                    lower = layer == 1 ? 0.0 : cutoff_radii[layer - 1]
                    upper = cutoff_radii[layer]
                    if lower < distance <= upper
                        push!(nth_neighbors[layer], j)
                        break
                    end
                end
            end
        end
        neighbors[i] = nth_neighbors
    end
    return neighbors
end

function get_ontop_sites(positions)
    [(p[1], p[2]) for p in positions]
end

function get_bridge_sites(positions::Vector{Vector{Float64}}, cutoff)
    cutoff2 = cutoff^2
    n = length(positions)

    sites = Vector{NTuple{2,Float64}}()

    for i in 1:n
        p1 = positions[i]

        for j in i+1:n
            p2 = positions[j]

            dx = p1[1]-p2[1]
            dy = p1[2]-p2[2]
            dz = p1[3]-p2[3]

            if dx*dx + dy*dy + dz*dz ≤ cutoff2
                push!(sites,
                      ((p1[1]+p2[1])/2,
                       (p1[2]+p2[2])/2))
            end
        end
    end

    return sites
end

function get_hollow_sites(positions, nn, tol)
    sites = Vector{NTuple{2,Float64}}()
    for p in positions
        right = nothing
        up    = nothing

        for q in positions
            dx = q[1]-p[1]
            dy = q[2]-p[2]
            dz = q[3]-p[3]

            d2 = dx*dx + dy*dy + dz*dz

            if abs(d2 - nn^2) ≤ tol
                if dx > tol && abs(dy) < tol
                    right = q
                elseif dy > tol && abs(dx) < tol
                    up = q
                end
            end
        end

        if right !== nothing && up !== nothing
            push!(sites,
                  ((p[1]+right[1])/2,
                   (p[2]+up[2])/2))
        end
    end

    return sites
end

function find_fcc_lattice_sites(slab, nn, tol)
    positions = get_positions(slab)
    ontop   = get_ontop_sites(positions)
    bridge  = get_bridge_sites(positions, nn)
    hollow  = get_hollow_sites(positions, nn, tol)

    return vcat(ontop, bridge, hollow)
end

function get_adsorbate_indicies(slab)
    adsorbate_indices = [i for i in 0:length(slab) - 1 if pyconvert(Int64, slab[i].tag) == 0]
    return adsorbate_indices
end

function get_adsorbate_positions(slab)
    adsorbate_indices = get_adsorbate_indicies(slab)
    adsorbate_positions = [slab[i].position for i in adsorbate_indices]
    return adsorbate_positions
end

function add_adsorbates!(slab, adsorbate_atoms, type_of_sites; height, coverage, nn, tol)
    positions = get_positions(slab)
    all_sites = []
    if "ontop" in type_of_sites
        ontop = get_ontop_sites(positions)
        all_sites = vcat(all_sites, ontop)
    end
    
    if "bridge" in type_of_sites
        bridge = get_bridge_sites(positions, nn)
        all_sites = vcat(all_sites, bridge)
    end
    
    if "hollow" in type_of_sites
        hollow = get_hollow_sites(positions, nn, tol)
        all_sites = vcat(all_sites, hollow)
    end
    n_ads = round(Int, coverage * length(all_sites))
    chosen = shuffle!(all_sites)[1:n_ads]
    adsorbate = adsorbate_atoms[1]
    for (x, y) in chosen
        ase.build.add_adsorbate(slab, adsorbate; height=height, position=(x, y))
    end

    adsorbate_indices = get_adsorbate_indicies(slab)
    return slab, all_sites, adsorbate_indices
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
                               components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                               adsorptions::Union{Vector{Int},Symbol}=:full)

    MLattice{C,TriangularLattice}(; lattice_constant::Float64=1.0,
                                  basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)],
                                  supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 2, 1),
                                  periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                  cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                  components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                  adsorptions::Union{Vector{Int},Symbol}=:full)

Constructs a square/triangular lattice with the specified parameters. The `components` and `adsorptions` arguments can be a vector of integers specifying
the indices of the occupied sites, or a symbol. If `components` is `:equal`, the lattice is divided into `C` equal components when possible, or 
nearest to equal components otherwise. If `adsorptions` is `:full`, all sites are classified as adsorption sites.

## Returns
- `MLattice{C,G}`: A square/triangular lattice object with `C` components.

"""
mutable struct MLattice{C,G} <: AbstractLattice
    lattice_vectors::Matrix{Float64}
    positions::Matrix{Float64}
    basis::Vector{Tuple{Float64, Float64, Float64}}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    periodicity::Tuple{Bool, Bool, Bool}
    cutoff_radii::Vector{Float64}
    components::Vector{Vector{Bool}}
    neighbors::Vector{Vector{Vector{Int}}}
    adsorptions::Vector{Bool}

    function MLattice{C,G}(
        lattice_vectors::Matrix{Float64},
        basis::Vector{Tuple{Float64, Float64, Float64}},
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        periodicity::Tuple{Bool, Bool, Bool},
        cutoff_radii::Vector{Float64},
        components::Vector{Vector{Bool}},
        adsorptions::Vector{Bool},
    ) where {C,G}

        num_components = length(components)

        if num_components != C
            throw(ArgumentError("For a $C-component system, got $num_components components!"))
        end

        positions = lattice_positions(lattice_vectors, basis, supercell_dimensions)

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii)
        
        return new{C,G}(lattice_vectors, positions, basis, supercell_dimensions, periodicity, cutoff_radii, components, neighbors, adsorptions)
    end
end


"""
    mutable struct AtomicLattice{C,G} <: AbstractLattice

A mutable struct representing an atomic lattice with adsorbates using ASE (Atomic Simulation Environment).

# Fields
- `lattice_atom::String`: The chemical symbol of the lattice substrate atom.
- `adsorbate_atoms::Vector{String}`: The chemical symbols of the adsorbate species.
- `coverage::Float64`: The fractional coverage of adsorbates on the surface.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `lattice_constant::Float64`: The lattice constant of the unit cell.
- `periodicity::Tuple{Bool, Bool, Bool}`: The periodic boundary conditions in each dimension.
- `lattice_positions::Matrix{Float64}`: The positions of the lattice points.
- `num_nearest_neighbors::Int64`: The number of nearest neighbors to consider.
- `neighbors::Vector{Vector{Vector{Int}}}`: The neighbor lists for each lattice point.
- `type_of_sites::Vector{String}`: The types of adsorption sites (e.g., "ontop", "bridge", "hollow").
- `ase_lattice::Py`: The ASE atoms object representing the complete system.
- `adsorbate_indices::Vector{Int64}`: The indices of adsorbate atoms in the ASE structure.
- `all_sites::Vector{Tuple{Float64, Float64}}`: The coordinates of all possible adsorption sites.

# Constructor
    AtomicLattice{C,G}(;
        lattice_atom::String,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String}=[""],
        coverage::Float64 = 0.5,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String}
    ) where {C,G}

Creates an `AtomicLattice` instance with the specified parameters. The constructor performs the following steps:
1. Validates that the number of adsorbate species matches the expected value `C`.
2. Constructs an FCC(100) slab using ASE with the specified lattice atom and dimensions.
3. Sets the periodic boundary conditions on the slab.
4. Adds adsorbates to the surface at the specified sites with the given coverage.
5. Computes the lattice positions and neighbor lists.

Throws an `ArgumentError` if the number of adsorbate species does not match `C`.

# Arguments
- `lattice_atom::String`: Chemical symbol for the substrate (e.g., "Pt", "Cu").
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: Size of supercell in (x, y, z) directions.
- `lattice_constant::Float64`: Lattice constant in Ångströms.
- `periodicity::Tuple{Bool, Bool, Bool}`: Periodic boundary conditions for each dimension.
- `adsorbate_atoms::Vector{String}`: Chemical symbols of adsorbates (default: `[""]`).
- `coverage::Float64`: Fractional surface coverage (default: `0.5`).
- `num_nearest_neighbors::Int64`: Number of nearest neighbors for neighbor list construction.
- `type_of_sites::Vector{String}`: Adsorption site types for each adsorbate species.

# Returns
- `AtomicLattice{C,G}`: An atomic lattice object with `C` adsorbate species and geometry type `G`.
"""

mutable struct AtomicLattice{C,G} <: AbstractLattice
    lattice_atom::String
    adsorbate_atoms::Vector{String}
    coverage::Float64
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    lattice_constant::Float64
    periodicity::Tuple{Bool, Bool, Bool}
    lattice_positions::Matrix{Float64}
    num_nearest_neighbors::Int64
    neighbors::Vector{Vector{Vector{Int}}}
    type_of_sites::Vector{String}
    ase_lattice::Py
    adsorbate_indices::Vector{Int64}
    all_sites::Vector{Tuple{Float64, Float64}}

    function AtomicLattice{C,G}(;
        lattice_atom::String,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String}=[""],
        coverage::Float64 = 0.5,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String}
    ) where {C,G}

        num_adsorbates = length(adsorbate_atoms)

        if num_adsorbates != C
            throw(ArgumentError("For a $C-adsorbate system, got $num_adsorbates adsorbates"))
        end

        slab = ase.build.fcc100(lattice_atom, supercell_dimensions, a=lattice_constant)
        slab.set_pbc(periodicity)
        
        if !isempty(adsorbate_atoms)
            ase_lattice, all_sites, adsorbate_indices = add_adsorbates!(slab, adsorbate_atoms, type_of_sites; height=1.0, coverage=coverage, nn=2.791, tol=0.1)
        end

        lattice_vectors = [lattice_constant 0 0; 0 lattice_constant 0; 0 0 1]
        lattice_positions = get_lattice_positions(lattice_vectors, supercell_dimensions)
        cutoff_radii = find_n_cutoff_radii(lattice_positions, num_nearest_neighbors)
        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, lattice_positions, periodicity, cutoff_radii)
        return new{C,G}(lattice_atom, adsorbate_atoms, coverage, supercell_dimensions, lattice_constant, periodicity, lattice_positions, num_nearest_neighbors, neighbors, type_of_sites, ase_lattice, adsorbate_indices, all_sites)
    end
end


"""
    split_into_subarrays(arr::AbstractVector, N::Int)

Split an array into `N` subarrays of approximately equal size.

# Arguments
- `arr::AbstractVector`: The array to split.
- `N::Int`: The number of subarrays to create.

# Returns
- `subarrays::Vector{Vector{eltype(arr)}}`: A vector of subarrays.

"""
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

"""
    mlattice_setup(C::Int, 
                     basis::Vector{Tuple{Float64, Float64, Float64}},
                     supercell_dimensions::Tuple{Int64, Int64, Int64},
                     components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol},
                     adsorptions::Union{Vector{Int}, Vector{Bool}, Symbol})

Setup the components and adsorptions for a lattice.

# Arguments
- `C::Int`: The number of components.
- `basis::Vector{Tuple{Float64, Float64, Float64}}`: The basis of the lattice.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}`: The components of the lattice.
- `adsorptions::Union{Vector{Int}, Vector{Bool}, Symbol}`: The adsorption sites on the lattice.

# Returns
- `lattice_comp::Vector{Vector{Bool}}`: The components of the lattice.
- `lattice_adsorptions::Vector{Bool}`: The adsorption sites on the lattice.

"""
function mlattice_setup(C::Int, 
                        basis::Vector{Tuple{Float64,Float64,Float64}},
                        supercell_dimensions::Tuple{Int64,Int64,Int64},
                        components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol},
                        adsorptions::Union{Vector{Int},Vector{Bool},Symbol})
    dim = prod(supercell_dimensions) * length(basis)
    lattice_adsorptions = zeros(Bool, dim)

    if adsorptions == :full
        lattice_adsorptions = [true for i in 1:dim]
    elseif adsorptions == :none
        lattice_adsorptions = [false for i in 1:dim]
    elseif adsorptions isa Vector{Int}
        for i in adsorptions
            lattice_adsorptions[i] = true
        end
    elseif adsorptions isa Vector{Bool}
        lattice_adsorptions = adsorptions
    else
        throw(ArgumentError("Adsorptions must be a vector of integers/booleans, or a supported symbol!"))
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
    elseif components isa Vector{Vector{Int}}
        lattice_comp = Vector{Vector{Bool}}(undef, C)
        for i in 1:C
            lattice_comp[i] = [false for i in 1:dim]
            for j in components[i]
                lattice_comp[i][j] = true
            end
        end
    elseif components isa Vector{Vector{Bool}}
        lattice_comp = components
    else
        throw(ArgumentError("components must be a vector of integers/booleans, or a supported symbol!"))
    end

    return lattice_comp, lattice_adsorptions
end

function MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
                                    basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
                                    supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
                                    periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                    cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                    components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                    adsorptions::Union{Vector{Int},Symbol}=:full,
                                ) where C

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 lattice_constant 0.0; 0.0 0.0 1.0]
    lattice_comp, lattice_adsorptions = mlattice_setup(C, basis, supercell_dimensions, components, adsorptions)

    return MLattice{C,SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, cutoff_radii, lattice_comp, lattice_adsorptions)
end

function MLattice{C,TriangularLattice}(; lattice_constant::Float64=1.0,
                                        basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)],
                                        supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 2, 1),
                                        periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                        cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                        components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                        adsorptions::Union{Vector{Int},Symbol}=:full,
                                    ) where C

    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 sqrt(3)*lattice_constant 0.0; 0.0 0.0 1.0]
    lattice_comp, lattice_adsorptions = mlattice_setup(C, basis, supercell_dimensions, components, adsorptions)

    return MLattice{C,TriangularLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, cutoff_radii, lattice_comp, lattice_adsorptions)
    
end





const SLattice{G} = MLattice{1,G} # alias for single-component lattices

const GLattice{C} = MLattice{C,GenericLattice} # alias for generic lattices

num_lattice_components(lattice::MLattice{C,G}) where {C,G} = C

"""
    num_sites(lattice::AbstractLattice)

Returns the total number of sites in a lattice given a `AbstractLattice` object. Returns the total number of sites.
"""
function num_sites(lattice::AbstractLattice)
    return prod(lattice.supercell_dimensions) * length(lattice.basis)
end

"""
    occupied_site_count(MLattice::MLattice{C})

Returns the number of occupied sites in each component of a lattice in an array.
"""
function occupied_site_count(MLattice::MLattice{C}) where C
    occupancy = Array{Int}(undef, C)
    for i in eachindex(MLattice.components)
        occupancy[i] = sum(MLattice.components[i])
    end
    return occupancy
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
mutable struct LatticeWalker{C} <: AbstractWalker
    configuration::AbstractLattice
    energy::typeof(0.0u"eV")
    iter::Int64
    function LatticeWalker(configuration::AbstractLattice; energy=0.0u"eV", iter=0)
        num_comp = num_lattice_components(configuration)
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