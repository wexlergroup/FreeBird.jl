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

get_positions(slab) = [pyconvert(Vector{Float64}, a.position) for a in slab]

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
    num_supercell_sites = prod(supercell_dimensions)

    a1, a2, a3 = [lattice_vectors[:, i] for i in 1:3]

    positions = Matrix{Float64}(undef, num_supercell_sites, 3)

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

"""
    compute_neighbors_banded(supercell_lattice_vectors::Matrix{Float64},
                             positions::Matrix{Float64},
                             periodicity::Tuple{Bool, Bool, Bool},
                             cutoff_radii::Vector{Float64})

Compute neighbour shells for an `AtomicLattice`, assigning each pair to the shell
whose *band* contains it: shell `k` holds `cutoff_radii[k-1] < d <= cutoff_radii[k]`
(with a lower bound of `0.0` for the first shell).

Named apart from `compute_neighbors` deliberately. It was introduced with
the same name **and the same argument types** as `compute_neighbors`, which on
Julia >= 1.12 is not a shadowing subtlety but a hard failure: precompiling the
package aborted with

    WARNING: Method definition compute_neighbors(...) in module AbstractWalkers
             at lattice_walkers.jl:20 overwritten at lattice_walkers.jl:170.
    ERROR: Method overwriting is not permitted during Module precompilation.

so `using FreeBird` could not load this branch at all.

How it differs from `compute_neighbors`, so the two can be reconciled later on
evidence rather than guesswork:

  1. Shell assignment. `compute_neighbors` takes the *first* shell whose cutoff
     the distance clears (`d <= cutoff_radii[i]`, then `break`). For cutoff radii
     in increasing order — which is what `find_n_cutoff_radii` produces — that is
     the same partition as the bands here, so on every current call site the two
     agree. They diverge only for unsorted `cutoff_radii`.
  2. Singular cells. This version builds a 2x2 reciprocal matrix when
     `!periodicity[3] && all(a3 .== 0)`, where `inv([a1 a2 a3])` would throw.
     Note the `AtomicLattice` constructor passes `a3 = [0, 0, 1] * nz`, so that
     branch is not reached from there today.
  3. Minimum image. This version applies the convention only when
     `periodicity[1] || periodicity[2]`, so a z-only-periodic cell would silently
     skip it; `compute_neighbors` always applies it (a no-op round trip when
     nothing is periodic). `compute_neighbors` is the more correct of the two here.

Test coverage pins claim 1 in `test/test-AbstractWalkers.jl`. If that equivalence
holds, the two functions should be collapsed into one — see MERGE_PLAN W10.
"""
function compute_neighbors_banded(supercell_lattice_vectors::Matrix{Float64}, 
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

"""
    add_adsorbates!(slab, adsorbate_atoms, type_of_sites; height, coverage, nn, tol)

Build the adsorption-site list for `slab` and decorate it to the requested
`coverage`. Returns `(slab, all_sites, occupations)`.

`occupations` is a `Vector{Bool}` over `all_sites` and is the ground truth for
which sites are filled. Two things about the previous version made that
impossible to express:

  * it called `shuffle!(all_sites)` **in place**, so the site list came back in
    random order with the occupied sites at the front. Occupancy was encoded as
    "the first `n_ads` entries of a list whose order is meaningless", which
    nothing downstream could read without knowing that;
  * it returned `adsorbate_indices` — indices into the *ASE frame* — which bear
    no relation to positions in `all_sites`. There was no way to ask whether
    site `i` was occupied.

The shuffle now happens in a permutation of the *indices*, `all_sites` keeps its
geometric order, and occupancy is a mask over it.
"""
function add_adsorbates!(slab, adsorbate_atoms, type_of_sites; height, coverage, nn, tol)
    positions = get_positions(slab)
    all_sites = Tuple{Float64,Float64}[]
    if "ontop" in type_of_sites
        append!(all_sites, get_ontop_sites(positions))
    end

    if "bridge" in type_of_sites
        append!(all_sites, get_bridge_sites(positions, nn))
    end

    if "hollow" in type_of_sites
        append!(all_sites, get_hollow_sites(positions, nn, tol))
    end

    isempty(all_sites) && throw(ArgumentError(
        "no adsorption sites produced for type_of_sites = $type_of_sites; " *
        "expected some of \"ontop\", \"bridge\", \"hollow\""))

    n_ads = round(Int, coverage * length(all_sites))
    occupations = fill(false, length(all_sites))
    occupations[randperm(length(all_sites))[1:n_ads]] .= true

    adsorbate = adsorbate_atoms[1]
    for i in findall(occupations)
        x, y = all_sites[i]
        ase.build.add_adsorbate(slab, adsorbate; height=height, position=(x, y))
    end

    return slab, all_sites, occupations
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
- `all_sites::Vector{Tuple{Float64, Float64}}`: Coordinates of every adsorption site, in geometric order.
- `occupations::Vector{Bool}`: **The ground truth.** Indexed over `all_sites`: site `i` is filled iff `occupations[i]`.
- `adsorbate_height::Float64`: Height at which adsorbates are placed. Stored rather than assumed, so the constructor and `sync_ase_lattice!` cannot disagree about it.
- `ase_dirty::Bool`: Whether `ase_lattice` is stale with respect to `occupations`.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `lattice_constant::Float64`: The lattice constant of the unit cell.
- `periodicity::Tuple{Bool, Bool, Bool}`: The periodic boundary conditions in each dimension.
- `lattice_positions::Matrix{Float64}`: The positions of the lattice points.
- `num_nearest_neighbors::Int64`: The number of nearest neighbors to consider.
- `neighbors::Vector{Vector{Vector{Int}}}`: The neighbor lists for each lattice point.
- `type_of_sites::Vector{String}`: The types of adsorption sites (e.g., "ontop", "bridge", "hollow").
- `ase_lattice::Py`: The ASE atoms object. A **derived cache** of `occupations`, not a second source of truth — see `sync_ase_lattice!`.

# Constructor
    AtomicLattice{C,G}(;
        lattice_atom::String,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String}=[""],
        coverage::Float64 = 0.5,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String},
        adsorbate_height::Float64 = 1.0
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
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    lattice_constant::Float64
    periodicity::Tuple{Bool, Bool, Bool}
    lattice_positions::Matrix{Float64}
    num_nearest_neighbors::Int64
    neighbors::Vector{Vector{Vector{Int}}}
    type_of_sites::Vector{String}
    all_sites::Vector{Tuple{Float64, Float64}}
    # ── ground truth ────────────────────────────────────────────────────────
    occupations::Vector{Bool}
    adsorbate_height::Float64
    # ── derived cache of the above; see sync_ase_lattice! ───────────────────
    ase_lattice::Py
    ase_dirty::Bool

    function AtomicLattice{C,G}(;
        lattice_atom::String,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String}=[""],
        coverage::Float64 = 0.5,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String},
        adsorbate_height::Float64 = 1.0
    ) where {C,G}

        num_adsorbates = length(adsorbate_atoms)

        if num_adsorbates != C
            throw(ArgumentError("For a $C-adsorbate system, got $num_adsorbates adsorbates"))
        end

        slab = ase.build.fcc100(lattice_atom, supercell_dimensions, a=lattice_constant)
        slab.set_pbc(periodicity)
        
        ase_lattice, all_sites, occupations = add_adsorbates!(
            slab, adsorbate_atoms, type_of_sites;
            height=adsorbate_height, coverage=coverage,
            nn=lattice_constant / sqrt(2), tol=0.1)

        lattice_vectors = [lattice_constant 0 0; 0 lattice_constant 0; 0 0 1]
        lattice_positions = get_lattice_positions(lattice_vectors, supercell_dimensions)
        cutoff_radii = find_n_cutoff_radii(lattice_positions, num_nearest_neighbors)
        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors_banded(supercell_lattice_vectors, lattice_positions, periodicity, cutoff_radii)
        return new{C,G}(lattice_atom, adsorbate_atoms, supercell_dimensions,
                        lattice_constant, periodicity, lattice_positions,
                        num_nearest_neighbors, neighbors, type_of_sites,
                        all_sites, occupations, adsorbate_height,
                        ase_lattice, false)
    end
end

"""
    coverage(lattice::AtomicLattice)

Fractional coverage, derived from `occupations`.

This was a stored field. It is computed now because a stored coverage is a
second copy of what `occupations` already says, and the two can disagree — the
whole point of W11 is that this type has one place where occupancy lives.
"""
coverage(lattice::AtomicLattice) = sum(lattice.occupations) / length(lattice.occupations)

"""
    nn_distance(lattice::AtomicLattice)

Surface nearest-neighbour distance, `a / sqrt(2)` for an fcc(100) termination.

The constructor previously passed a hardcoded `nn = 2.791` to `add_adsorbates!`.
That is `3.947 / sqrt(2)` — the value for palladium — so the site-finding
geometry was silently correct for exactly one `lattice_constant` and wrong for
every other, with no error, just a different (or empty) set of adsorption sites.
"""
nn_distance(lattice::AtomicLattice) = lattice.lattice_constant / sqrt(2)

"""
    sync_ase_lattice!(lattice::AtomicLattice)

Rebuild `ase_lattice`'s adsorbates from `occupations`, and clear `ase_dirty`.

`ase_lattice` is a cache. Occupancy moves update `occupations` and set
`ase_dirty`; the ASE frame is only made to agree when something actually needs
to look at it — an energy evaluation through a Python calculator, or writing a
trajectory. Doing it eagerly on every move would put a Python round trip in the
innermost Monte Carlo loop, which at ~16k proposals per walk is where all the
time would go.

Which means: **anything that reads `ase_lattice` must call this first.** That is
the one rule this design imposes, and it is the reason the dirty flag is a field
rather than a convention.

Adsorbates are identified by ASE tag 0, which is what `ase.build.add_adsorbate`
assigns and what `get_adsorbate_indicies` already relies on; the substrate keeps
the layer tags `fcc100` gave it. Deletion goes in reverse index order because
removing an atom renumbers everything after it.
"""
function sync_ase_lattice!(lattice::AtomicLattice)
    lattice.ase_dirty || return lattice

    slab = lattice.ase_lattice
    tags = pyconvert(Vector{Int}, slab.get_tags())
    ads = findall(==(0), tags)
    if !isempty(ads)
        slab.__delitem__(pylist([i - 1 for i in reverse(ads)]))
    end

    adsorbate = lattice.adsorbate_atoms[1]
    for i in findall(lattice.occupations)
        x, y = lattice.all_sites[i]
        ase.build.add_adsorbate(slab, adsorbate;
                                height=lattice.adsorbate_height, position=(x, y))
    end

    lattice.ase_dirty = false
    return lattice
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
    num_lattice_components(lattice::AtomicLattice{C,G}) where {C,G}

Number of adsorbate species on an `AtomicLattice`, i.e. its first type parameter.

Small, but it is what makes `AtomicLattice` usable as a walker configuration at
all: `LatticeWalker`'s inner constructor calls `num_lattice_components` to fix
its own type parameter, so without a method here every
`LatticeWalker(::AtomicLattice)` is a `MethodError`.
"""
num_lattice_components(lattice::AtomicLattice{C,G}) where {C,G} = C

"""
    num_sites(lattice::AbstractLattice)

Returns the total number of sites in a lattice given a `AbstractLattice` object. Returns the total number of sites.
"""
function num_sites(lattice::AbstractLattice)
    return prod(lattice.supercell_dimensions) * length(lattice.basis)
end

"""
    num_sites(lattice::AtomicLattice)

Number of adsorption sites on an `AtomicLattice`: the length of `all_sites`.

The generic `AbstractLattice` method above cannot serve here. It computes
`prod(supercell_dimensions) * length(basis)` — the substrate grid — and
`AtomicLattice` has no `basis` field at all, so the generic method is a
`FieldError` rather than a wrong answer. The two counts are also genuinely
different: `all_sites` is a union of ontop, bridge and hollow positions selected
by `type_of_sites`, and is not the substrate grid.
"""
num_sites(lattice::AtomicLattice) = length(lattice.all_sites)

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