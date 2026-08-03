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

"""
    _minimum_image_distance(supercell_lattice_vectors, reciprocal_lattice_vectors,
                            periodicity, pos_i, pos_j)

Minimum-image distance between two Cartesian positions under the given
supercell and periodicity — the single distance kernel shared by
`compute_neighbors` and `enumerate_motif_embeddings`, so pair shells and
cluster embeddings follow the identical torus convention.
"""
@inline function _minimum_image_distance(supercell_lattice_vectors::AbstractMatrix{Float64},
                                         reciprocal_lattice_vectors::AbstractMatrix{Float64},
                                         periodicity::Tuple{Bool,Bool,Bool},
                                         pos_i, pos_j)
    dr = [pos_j[1] - pos_i[1], pos_j[2] - pos_i[2], pos_j[3] - pos_i[3]]
    fractional_dr = reciprocal_lattice_vectors * dr
    for k in 1:3
        if periodicity[k]
            fractional_dr[k] -= round(fractional_dr[k])
        end
    end
    return norm(supercell_lattice_vectors * fractional_dr)
end

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

        for j in 1:num_atoms
            if i != j
                distance = _minimum_image_distance(
                    supercell_lattice_vectors, reciprocal_lattice_vectors,
                    periodicity, view(positions, i, :), view(positions, j, :))

                for k in 1:layers_of_neighbors
                    if distance <= cutoff_radii[k]
                        push!(nth_neighbors[k], j)
                        break
                    end
                end
            end
        end

        neighbors[i] = nth_neighbors
    end

    for k in 1:layers_of_neighbors
        if all(isempty(nbrs[k]) for nbrs in neighbors)
            @warn "Neighbor shell $k (cutoff radius $(cutoff_radii[k])) is empty on every " *
                  "site; a Hamiltonian coupling for this shell silently contributes zero. " *
                  "Check cutoff_radii against the actual neighbor distances — e.g. on a " *
                  "triangular lattice the second-neighbor distance is √3 ≈ 1.732 lattice " *
                  "constants, beyond the square-lattice cutoff of 1.5 — or whether the " *
                  "supercell is too small to contain this shell."
        end
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
                                  cutoff_radii::Vector{Float64}=[1.1, 1.8],
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
                                        # The triangular second-neighbor distance is √3 ≈ 1.732, so the
                                        # square-lattice default [1.1, 1.5] would leave shell 2 empty.
                                        cutoff_radii::Vector{Float64}=[1.1, 1.8],
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
    order_parameter_c2x2(lattice::MLattice{1,SquareLattice}) -> Float64

Sublattice order parameter for c(2×2) checkerboard ordering on a
single-component square lattice:

    Ψ = |Σ_{occupied sites} (−1)^(i+j)| / M,

where `(i, j)` are the integer lattice coordinates of each site and `M` is
the number of sites. Equivalent to the four-sublattice
`Ψ_c(2×2) = |N_a + N_d − N_b − N_c| / M` of Zhang, Blum & Reuter
[PRB **75**, 235406 (2007), Eq. (8)]: sublattices (a, d) share one
checkerboard parity and (b, c) the other. A perfect c(2×2) arrangement at
half filling gives `1/2`; the empty and the full lattice give `0`.

Ψ is not a function of `(E, N)` — degenerate energy levels contain ordered
and disordered configurations alike — so it must be evaluated on
configurations (e.g. per culled walker, via the `observables` keyword of the
nested-sampling loops) rather than reconstructed from an energy ledger.

Requires a single-site basis, a strictly two-dimensional supercell
(`supercell_dimensions[3] == 1`), and even in-plane dimensions (the two
checkerboard sublattices must tile the periodic cell evenly); violations
throw an `ArgumentError`.
"""
function order_parameter_c2x2(lattice::MLattice{1,SquareLattice})
    if length(lattice.basis) != 1
        throw(ArgumentError("order_parameter_c2x2 requires a single-site " *
            "basis, got $(length(lattice.basis)) basis sites"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("order_parameter_c2x2 requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    if isodd(d1) || isodd(d2)
        throw(ArgumentError("order_parameter_c2x2 requires even in-plane " *
            "supercell dimensions so the checkerboard sublattices tile the " *
            "periodic cell evenly, got ($d1, $d2)"))
    end
    occ = lattice.components[1]
    acc = 0
    for s in eachindex(occ)
        if occ[s]
            # lattice_positions ordering: basis innermost, dimension 1 fastest
            i0 = (s - 1) % d1
            j0 = ((s - 1) ÷ d1) % d2
            acc += iseven(i0 + j0) ? 1 : -1
        end
    end
    return abs(acc) / (d1 * d2)
end

"""
    order_parameter_sqrt3(lattice::MLattice{1,TriangularLattice}) -> Float64

Three-sublattice order parameter for (√3×√3)R30° ordering on a
single-component triangular lattice:

    Ψ = |Σ_{occupied sites} ω^{c(s)}| / M,   ω = e^{2πi/3},

where `c(s) ∈ {0,1,2}` is the site's sublattice label under the standard
tripartition of the triangular lattice and `M` is the number of sites. This
is the modulus of the complex three-state Potts order parameter: the Z₃
phase distinguishing the three degenerate ordered states is divided out, so
all three give the same value. A perfect √3×√3 arrangement at coverage 1/3
gives `1/3`; the empty and the full lattice give `0` (1 + ω + ω² = 0).

Ψ is not a function of `(E, N)`: degenerate energy levels contain ordered
and disordered configurations alike, so it must be evaluated on
configurations (e.g. per culled walker, via the `observables` keyword of
the nested-sampling loops). Higher moments for Binder-cumulant analysis are
composed caller-side, with no library change:

    observables = [:psi  => order_parameter_sqrt3,
                   :psi2 => cfg -> order_parameter_sqrt3(cfg)^2,
                   :psi4 => cfg -> order_parameter_sqrt3(cfg)^4]

after which ⟨Ψ⟩, ⟨Ψ²⟩, ⟨Ψ⁴⟩ come from `observable_cols=[:psi, :psi2, :psi4]`
in the grand-canonical stats functions.

Requires the standard two-site centered-rectangular triangular basis
`[(0, 0, 0), (a/2, √3·a/2, 0)]` consistent with the lattice vectors, a
strictly two-dimensional supercell (`supercell_dimensions[3] == 1`), and
`supercell_dimensions[1]` divisible by 3 (the tripartition closes on the
periodic cell iff the a₁ circumference is a multiple of 3; the a₂ dimension
is unconstrained); violations throw an `ArgumentError`. Note the shipped
default `supercell_dimensions = (4, 2, 1)` is *not* commensurate.
"""
function order_parameter_sqrt3(lattice::MLattice{1,TriangularLattice})
    if length(lattice.basis) != 2
        throw(ArgumentError("order_parameter_sqrt3 requires the two-site " *
            "centered-rectangular triangular basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    ax, ay = lattice.lattice_vectors[1, 1], lattice.lattice_vectors[2, 2]
    b1, b2 = lattice.basis
    if !(isapprox(b1[1], 0.0, atol=1e-9) && isapprox(b1[2], 0.0, atol=1e-9) &&
         isapprox(b2[1], ax / 2, atol=1e-9) && isapprox(b2[2], ay / 2, atol=1e-9))
        throw(ArgumentError("order_parameter_sqrt3 requires the standard " *
            "triangular basis [(0, 0, 0), (a/2, √3·a/2, 0)] consistent with " *
            "the lattice vectors, got $(lattice.basis)"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("order_parameter_sqrt3 requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    if d1 % 3 != 0
        throw(ArgumentError("order_parameter_sqrt3 requires " *
            "supercell_dimensions[1] divisible by 3 so the three √3×√3 " *
            "sublattices close on the periodic cell, got $d1"))
    end
    occ = lattice.components[1]
    # Sublattice occupations n[c+1], c ∈ {0,1,2}. lattice_positions ordering:
    # basis innermost, dimension 1 next. In triangular integer coordinates
    # (m, n) with r = m·t₁ + n·t₂, a valid tripartition is c = (m − n) mod 3;
    # the basis-0 site of cell column ci sits at (ci − cj, 2cj), giving
    # c = ci mod 3, and the basis-1 site adds one t₂ step, giving
    # c = (ci + 2) mod 3. Wrapping the a₁ circumference shifts c by d1 mod 3,
    # hence the commensurability guard above.
    n = zeros(Int, 3)
    for s in eachindex(occ)
        occ[s] || continue
        b = (s - 1) % 2
        ci = ((s - 1) ÷ 2) % d1
        c = b == 0 ? ci % 3 : (ci + 2) % 3
        n[c+1] += 1
    end
    M = 2 * d1 * d2
    # |n₀ + n₁·ω + n₂·ω²|² = n₀² + n₁² + n₂² − n₀n₁ − n₁n₂ − n₂n₀, exact in
    # integers; one final square root, no complex arithmetic.
    s2 = n[1]^2 + n[2]^2 + n[3]^2 - n[1] * n[2] - n[2] * n[3] - n[3] * n[1]
    return sqrt(Float64(s2)) / M
end

"""
    motif_distances(coords::AbstractVector{<:Tuple}) -> Vector{Float64}

Sorted multiset of the pairwise Euclidean distances of a coordinate
template, for transcribing a cluster-interaction figure straight from a
paper's geometry (e.g. the trio `[(0, 0), (1, 0), (0, 1)]` gives
`[1, 1, √2]`). Two-component tuples are treated as in-plane coordinates
with `z = 0`. The result is the `distances` argument of
[`enumerate_motif_embeddings`](@ref).
"""
function motif_distances(coords::AbstractVector{<:Tuple})
    n = length(coords)
    n >= 2 || throw(ArgumentError("a motif needs at least two sites"))
    to3(t) = (Float64(t[1]), Float64(t[2]), length(t) >= 3 ? Float64(t[3]) : 0.0)
    pts = [to3(t) for t in coords]
    ds = Float64[]
    for a in 1:n, b in (a+1):n
        push!(ds, sqrt((pts[a][1] - pts[b][1])^2 +
                       (pts[a][2] - pts[b][2])^2 +
                       (pts[a][3] - pts[b][3])^2))
    end
    return sort(ds)
end

"""
    enumerate_motif_embeddings(lattice::MLattice, distances::AbstractVector{<:Real};
                               tol::Float64=1e-6,
                               expected_count::Union{Int,Nothing}=nothing)
        -> Vector{NTuple{K,Int}}

Enumerate every embedding of a cluster motif on a periodic lattice. The
motif is declared by the sorted multiset of its pairwise minimum-image
distances (`K` is inferred from the multiset length: 1 → pair, 3 → trio,
6 → quattro, 10 → quinto); [`motif_distances`](@ref) builds the multiset
from a coordinate template. Distances are compared with absolute tolerance
`tol`, using the same minimum-image kernel as `compute_neighbors`, so
embeddings follow the identical torus convention as the pair shells: one
entry per unordered site set whose distances match, in canonical strictly
increasing order — the form [`ClusterInteraction`](@ref) requires.

Diagnostics:
- The total embedding count and the count per site are logged (`@info`).
- A warning is emitted when the per-site embedding membership is not
  uniform: on a site-transitive lattice, nonuniformity indicates distance
  aliasing or an unsuitable `tol`.
- A warning is emitted when any periodic circumference does not exceed
  `K · maximum(distances)`: such a cell is not a faithful quotient, and
  winding (wrap-around) embeddings are counted under the torus convention
  (e.g. the 18-site triangular cell carries 42 nearest-neighbor-triangle
  embeddings: 36 faces plus 6 winding three-cycles).
- When `expected_count` is given (e.g. a hand-derived per-cell
  multiplicity), a mismatch throws an `ArgumentError` — the recommended
  guard against silent transcription errors.

Note: for a pair signature (`K = 2`) this reproduces the sites of a
neighbor shell as unordered pairs — useful as a counting diagnostic; pair
couplings themselves belong in `GenericLatticeHamiltonian`.

**Homometry caveat**: for `K ≥ 4`, non-congruent figures can share a
distance multiset (homometric figures), and this method then enumerates
the embeddings of every such figure together — it warns about this. Pass
the coordinate template instead (the `coords` method) to enumerate only
embeddings whose full distance *matrix* matches the template under some
site permutation, which excludes homometric aliases while preserving the
torus counting convention. For `K ≤ 3` the multiset determines the figure
and the two methods agree.
"""
function enumerate_motif_embeddings(lattice::MLattice, distances::AbstractVector{<:Real};
                                    tol::Float64=1e-6,
                                    expected_count::Union{Int,Nothing}=nothing)
    npairs = length(distances)
    K = round(Int, (1 + sqrt(1 + 8.0 * npairs)) / 2)
    if K * (K - 1) ÷ 2 != npairs || !(2 <= K <= 5)
        throw(ArgumentError(
            "distances must be the full pairwise multiset of a K-site motif " *
            "— 1 (K = 2), 3 (K = 3), 6 (K = 4), or 10 (K = 5) entries — got $npairs"))
    end
    if K >= 4
        @warn "a distance multiset does not determine a figure for K ≥ 4 " *
              "(homometric figures share multisets); embeddings of every " *
              "figure with this multiset are enumerated together. Pass the " *
              "coordinate template to enumerate_motif_embeddings to select " *
              "the congruent embeddings only."
    end
    return _enumerate_motif_core(lattice, collect(Float64, distances), K,
                                 nothing, tol, expected_count)
end

"""
    enumerate_motif_embeddings(lattice::MLattice, coords::AbstractVector{<:Tuple};
                               tol=1e-6, expected_count=nothing)

Template method: declare the motif by its site coordinates (as accepted by
[`motif_distances`](@ref)) and enumerate only the embeddings whose full
minimum-image distance matrix matches the template's under some site
permutation. This is the recommended method for `K ≥ 4`, where a distance
multiset alone does not determine the figure (homometric figures).
"""
function enumerate_motif_embeddings(lattice::MLattice, coords::AbstractVector{<:Tuple};
                                    tol::Float64=1e-6,
                                    expected_count::Union{Int,Nothing}=nothing)
    K = length(coords)
    2 <= K <= 5 || throw(ArgumentError(
        "the motif template must have 2 to 5 sites, got $K"))
    to3(t) = (Float64(t[1]), Float64(t[2]), length(t) >= 3 ? Float64(t[3]) : 0.0)
    pts = [to3(t) for t in coords]
    T = zeros(K, K)
    for a in 1:K, b in 1:K
        T[a, b] = sqrt((pts[a][1] - pts[b][1])^2 +
                       (pts[a][2] - pts[b][2])^2 +
                       (pts[a][3] - pts[b][3])^2)
    end
    sig = sort([T[a, b] for a in 1:K for b in (a+1):K])
    return _enumerate_motif_core(lattice, sig, K, T, tol, expected_count)
end

# Does some permutation of `sites` match the template distance matrix `T`
# entrywise within tol? Backtracking assignment with early pruning; K ≤ 5,
# so at most 120 permutations are ever considered.
function _matches_template(sites::Vector{Int}, dist, T::Matrix{Float64}, tol::Float64)
    K = length(sites)
    perm = zeros(Int, K)
    used = falses(K)
    function assign(a)
        a > K && return true
        for b in 1:K
            used[b] && continue
            ok = true
            for a2 in 1:(a-1)
                if abs(dist(sites[perm[a2]], sites[b]) - T[a2, a]) > tol
                    ok = false
                    break
                end
            end
            ok || continue
            used[b] = true
            perm[a] = b
            assign(a + 1) && return true
            used[b] = false
        end
        return false
    end
    return assign(1)
end

function _enumerate_motif_core(lattice::MLattice, distances::Vector{Float64}, K::Int,
                               template::Union{Nothing,Matrix{Float64}},
                               tol::Float64, expected_count::Union{Int,Nothing})
    npairs = length(distances)
    sig = sort(distances)
    all(>(0.0), sig) || throw(ArgumentError("motif distances must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    M = num_sites(lattice)
    dims = lattice.supercell_dimensions
    scv = lattice.lattice_vectors * Diagonal([dims[1], dims[2], dims[3]])
    rec = inv(scv)
    pos = lattice.positions
    per = lattice.periodicity

    # Wrap-around guard: the shortest periodic translation (over small
    # integer combinations of the periodic supercell vectors, which covers
    # reasonably reduced cells; extreme shear beyond ±1 combinations is not
    # detected) must exceed K·d_max, else winding embeddings are admitted
    # and counts follow the torus convention
    d_max = sig[end]
    C_min = Inf
    for n1 in -1:1, n2 in -1:1, n3 in -1:1
        (n1 == 0 && n2 == 0 && n3 == 0) && continue
        (!per[1] && n1 != 0) && continue
        (!per[2] && n2 != 0) && continue
        (!per[3] && n3 != 0) && continue
        C_min = min(C_min, norm(n1 * scv[:, 1] + n2 * scv[:, 2] + n3 * scv[:, 3]))
    end
    if isfinite(C_min) && C_min <= K * d_max + tol
        @warn "shortest periodic translation $(round(C_min, digits=4)) does " *
              "not exceed K·d_max = $(K * d_max); the cell is not a faithful " *
              "quotient and embeddings follow the torus (minimum-image) " *
              "convention"
    end

    uniq = Float64[]
    for d in sig
        any(u -> abs(u - d) <= tol, uniq) || push!(uniq, d)
    end

    dist(i, j) = _minimum_image_distance(scv, rec, per,
                                         view(pos, i, :), view(pos, j, :))

    # Per-site candidates at any signature distance, restricted to j > i:
    # anchored strictly increasing enumeration makes each unordered set
    # appear exactly once, with the anchor as its minimum index. Pruning
    # measures 2·tol from the merged uniq representatives so it stays a
    # relaxation of the elementwise acceptance below (a representative can
    # sit up to tol away from the signature entry it absorbed).
    cand = [Int[] for _ in 1:M]
    for i in 1:M, j in (i+1):M
        d = dist(i, j)
        if any(u -> abs(u - d) <= 2 * tol, uniq)
            push!(cand[i], j)
        end
    end

    embeddings = Vector{NTuple{K,Int}}()
    partial = Int[]
    function extend!()
        if length(partial) == K
            ds = Float64[]
            for a in 1:K, b in (a+1):K
                push!(ds, dist(partial[a], partial[b]))
            end
            sort!(ds)
            if all(abs(ds[t] - sig[t]) <= tol for t in 1:npairs)
                if template === nothing || _matches_template(partial, dist, template, tol)
                    push!(embeddings, NTuple{K,Int}(partial))
                end
            end
            return nothing
        end
        for j in cand[partial[1]]
            j > partial[end] || continue
            ok = true
            for s in partial
                dj = dist(s, j)
                if !any(u -> abs(u - dj) <= 2 * tol, uniq)
                    ok = false
                    break
                end
            end
            ok || continue
            push!(partial, j)
            extend!()
            pop!(partial)
        end
        return nothing
    end
    for i in 1:M
        empty!(partial)
        push!(partial, i)
        extend!()
    end

    counts = zeros(Int, M)
    for e in embeddings, s in e
        counts[s] += 1
    end
    total = length(embeddings)
    @info "enumerate_motif_embeddings: $total embeddings of a $K-site motif " *
          "($(round(total / M, digits=4)) per site)"
    if !isempty(embeddings) && !all(==(counts[1]), counts)
        @warn "per-site embedding membership is not uniform " *
              "(min $(minimum(counts)), max $(maximum(counts))); on a " *
              "site-transitive lattice this indicates distance aliasing or " *
              "an unsuitable tol"
    end
    if expected_count !== nothing && total != expected_count
        throw(ArgumentError(
            "enumerated $total embeddings but expected_count = " *
            "$expected_count; check the motif signature, tol, and the cell " *
            "size (wrap-around)"))
    end
    return embeddings
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