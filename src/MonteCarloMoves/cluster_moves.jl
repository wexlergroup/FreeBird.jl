"""
    geometric_cluster_swap!(lattice::MLattice{C,SquareLattice}, p::Float64) where C

Perform a geometric cluster move on a square lattice by point inversion through a random pivot.

A random pivot site and seed site are chosen on the periodic lattice. A cluster of sites is
grown from the seed via BFS, adding first-shell neighbors with fixed probability `p`
(independent of the configuration). Each site in the cluster is paired with its reflection
through the pivot, and the occupation states of each pair are swapped across all components.

Because cluster growth does not depend on the occupancy field, the proposal is symmetric:
P(cluster | config_A) = P(cluster | config_B). This means no Metropolis correction is needed;
in nested sampling, one simply checks E_proposed < E_max.

The parameter `p` controls the average cluster size: `p ≈ 0` gives single-site moves,
`p → 1` gives lattice-spanning clusters.

# Arguments
- `lattice::MLattice{C,SquareLattice}`: The lattice to perform the cluster move on.
- `p::Float64`: Growth probability for BFS cluster construction (0 < p < 1).
- `record::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing`: When a vector is
  supplied, every applied (site, reflected_site) pair is appended to it;
  re-applying the recorded pairs (each pair swap is an involution) reverts
  the move exactly. The `nothing` default changes no behavior and no
  random-number draw.

# Returns
- `lattice::MLattice{C,SquareLattice}`: The mutated lattice after the cluster swap.

# Notes
- Preserves per-component particle counts exactly.
- The move is self-inverse for a fixed cluster: applying the same cluster swap twice
  restores the original configuration.
- Supports 2D (`supercell_dimensions[3] == 1`) and 3D lattices. In 3D, the reflection
  operates in the periodic xy-plane; the z-coordinate is preserved.
- Requires the single-site basis (the site-index reflection assumes single-basis,
  x-fastest ordering); a multi-site basis throws an `ArgumentError`.

# References
- Heringa & Blöte, Phys. Rev. E 57, 4976 (1998) — geometric cluster MC framework.
- Adaptation uses fixed growth probability (configuration-independent) for symmetric proposals.
"""
function geometric_cluster_swap!(lattice::MLattice{C,SquareLattice}, p::Float64;
                                 record::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing) where C
    if length(lattice.basis) != 1
        throw(ArgumentError("geometric_cluster_swap! on a square lattice " *
            "requires the single-site basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    Lx = lattice.supercell_dimensions[1]
    Ly = lattice.supercell_dimensions[2]
    Lz = lattice.supercell_dimensions[3]
    n_sites = num_sites(lattice)

    # Choose random pivot (grid coordinates, 0-indexed) in the periodic xy-plane
    pivot_gx = rand(0:Lx-1)
    pivot_gy = rand(0:Ly-1)

    # Choose random seed site (1-indexed)
    seed = rand(1:n_sites)

    # Build cluster via BFS with fixed growth probability
    reflect = s -> _reflect_site(s, pivot_gx, pivot_gy, Lx, Ly, Lz)
    cluster = _build_geometric_cluster(lattice, seed, reflect, p)

    # Swap occupation states for each applied (site, reflected_site) pair;
    # the recorded pairs revert the move when re-applied (involution)
    _apply_cluster_pairs!(lattice, cluster)
    if record !== nothing
        for pr in cluster
            pr[1] != pr[2] && push!(record, pr)
        end
    end

    return lattice
end

# --- Internal helpers (not exported) ---

"""
    _apply_cluster_pairs!(lattice::MLattice{C,G}, pairs::Vector{Tuple{Int,Int}}) where {C,G}

Apply the occupancy exchange for each (site, reflected_site) pair with distinct
indices, across all components. Each pair swap is an involution, so re-applying
the same pair list reverts a cluster move exactly; the copy-free kernels use it
as the revert path for rejected cluster proposals.
"""
function _apply_cluster_pairs!(lattice::MLattice{C,G},
                               pairs::Vector{Tuple{Int,Int}}) where {C,G}
    for (a, b) in pairs
        if a != b
            for comp in 1:C
                lattice.components[comp][a], lattice.components[comp][b] =
                    lattice.components[comp][b], lattice.components[comp][a]
            end
        end
    end
    return lattice
end

"""
    _site_to_grid(site::Int, Lx::Int, Ly::Int) -> Tuple{Int,Int,Int}

Convert a 1-indexed site index to 0-indexed 3D grid coordinates (gx, gy, gz)
for a single-basis square lattice. Site ordering: x varies fastest, then y, then z.
"""
@inline function _site_to_grid(site::Int, Lx::Int, Ly::Int)
    s = site - 1
    gz = s ÷ (Lx * Ly)
    rem_xy = s % (Lx * Ly)
    gy = rem_xy ÷ Lx
    gx = rem_xy % Lx
    return (gx, gy, gz)
end

"""
    _grid_to_site(gx::Int, gy::Int, gz::Int, Lx::Int, Ly::Int) -> Int

Convert 0-indexed 3D grid coordinates to a 1-indexed site index.
"""
@inline function _grid_to_site(gx::Int, gy::Int, gz::Int, Lx::Int, Ly::Int)
    return gz * Lx * Ly + gy * Lx + gx + 1
end

"""
    _reflect_site(site::Int, pivot_gx::Int, pivot_gy::Int, Lx::Int, Ly::Int, Lz::Int) -> Int

Compute the point inversion of `site` through the pivot at grid coordinates
`(pivot_gx, pivot_gy)` with periodic wrapping in x and y. The z-coordinate
is preserved (no reflection in the non-periodic z direction).
"""
@inline function _reflect_site(site::Int, pivot_gx::Int, pivot_gy::Int, Lx::Int, Ly::Int, Lz::Int)
    gx, gy, gz = _site_to_grid(site, Lx, Ly)
    rx = mod(2 * pivot_gx - gx, Lx)
    ry = mod(2 * pivot_gy - gy, Ly)
    return _grid_to_site(rx, ry, gz, Lx, Ly)
end

"""
    _build_geometric_cluster(lattice, seed, reflect, p) -> Vector{Tuple{Int,Int}}

Grow a geometric cluster from `seed` using BFS with fixed growth probability `p`,
pairing each visited site with its image under the `reflect` closure (a
self-inverse site-index map supplied by the geometry-specific caller).
Returns a vector of (site, reflected_site) pairs.
"""
function _build_geometric_cluster(lattice::MLattice,
                                  seed::Int,
                                  reflect::F,
                                  p::Float64) where {F}
    visited = Set{Int}()
    stack = Int[seed]
    cluster = Tuple{Int,Int}[]

    while !isempty(stack)
        s = pop!(stack)

        s in visited && continue

        push!(visited, s)
        r = reflect(s)
        push!(visited, r)
        push!(cluster, (s, r))

        # Grow to first-shell neighbors with probability p
        for n in lattice.neighbors[s][1]
            if n ∉ visited && rand() < p
                push!(stack, n)
            end
        end
    end

    return cluster
end

# --- Triangular lattice (two-site centered-rectangular basis) ---

"""
    _tri_site_to_halfgrid(site::Int, nx::Int) -> Tuple{Int,Int}

Half-grid coordinates (hx, hy) = (2·ci + b, 2·cj + b) of a site on the
two-site-basis triangular lattice, in units of a/2 and √3·a/2, where
b = basis index and (ci, cj) the 0-indexed cell. Under the standard site
ordering (basis innermost, dimension 1 next), the site set is exactly
{(hx, hy) : hx ≡ hy (mod 2)} on the (2·nx, 2·ny) torus — the
centered-rectangular condition.
"""
@inline function _tri_site_to_halfgrid(site::Int, nx::Int)
    s = site - 1
    b = s % 2
    cell = s ÷ 2
    ci = cell % nx
    cj = cell ÷ nx
    return (2 * ci + b, 2 * cj + b)
end

"""
    _tri_halfgrid_to_site(hx::Int, hy::Int, nx::Int) -> Int

Inverse of [`_tri_site_to_halfgrid`](@ref) for in-range half-grid
coordinates with hx ≡ hy (mod 2).
"""
@inline function _tri_halfgrid_to_site(hx::Int, hy::Int, nx::Int)
    b = hx & 1
    ci = (hx - b) >> 1
    cj = (hy - b) >> 1
    return b + 2 * (ci + nx * cj) + 1
end

"""
    _tri_reflect_site(site::Int, hpx::Int, hpy::Int, nx::Int, ny::Int) -> Int

Point inversion of `site` through the pivot whose doubled half-grid
position is (hpx, hpy) = h(s₁) + h(s₂) for two pivot sites s₁, s₂, with
periodic wrapping: σ(h) = (hpx, hpy) − h mod (2·nx, 2·ny). Both pivot
sites satisfy hx ≡ hy (mod 2), so hpx ≡ hpy (mod 2) and σ preserves the
parity difference — the reflected point is always a site. Site pivots
(s₁ = s₂, hpx even) preserve each site's basis index; midpoint pivots
with hpx odd exchange the two basis sublattices. σ is an involution:
σ(σ(h)) = h.
"""
@inline function _tri_reflect_site(site::Int, hpx::Int, hpy::Int, nx::Int, ny::Int)
    hx, hy = _tri_site_to_halfgrid(site, nx)
    rx = mod(hpx - hx, 2 * nx)
    ry = mod(hpy - hy, 2 * ny)
    return _tri_halfgrid_to_site(rx, ry, nx)
end

"""
    _tri_site_to_halfgrid(site::Int, nx::Int, ny::Int) -> Tuple{Int,Int,Int}

Layer-resolved half-grid coordinates (hx, hy, k) of a site on a stacked
two-site-basis triangular cell, with k the 0-based layer index. Under the
standard site ordering (basis innermost, dimension 1 next, dimension 3
outermost) the sites of layer k occupy the contiguous block
`(k B + 1):((k + 1) B)`, B = 2 nx ny, and the in-plane coordinates are those
of the two-argument [`_tri_site_to_halfgrid`](@ref) applied to the
block-local index. At nz = 1 the layer index is identically zero and the
in-plane coordinates are exactly those of the two-argument method.
"""
@inline function _tri_site_to_halfgrid(site::Int, nx::Int, ny::Int)
    k, s0 = divrem(site - 1, 2 * nx * ny)
    hx, hy = _tri_site_to_halfgrid(s0 + 1, nx)
    return (hx, hy, k)
end

"""
    _tri_halfgrid_to_site(hx::Int, hy::Int, k::Int, nx::Int, ny::Int) -> Int

Inverse of the three-argument [`_tri_site_to_halfgrid`](@ref) for
in-range half-grid coordinates with hx ≡ hy (mod 2) and layer index
0 ≤ k < nz.
"""
@inline function _tri_halfgrid_to_site(hx::Int, hy::Int, k::Int, nx::Int, ny::Int)
    return k * (2 * nx * ny) + _tri_halfgrid_to_site(hx, hy, nx)
end

"""
    _tri_reflect_site_3d(site::Int, hpx::Int, hpy::Int, kp::Int, nx::Int, ny::Int, nz::Int) -> Int

Point inversion of `site` through the midpoint of two pivot sites
(h₁, k₁) and (h₂, k₂) on a stacked two-site-basis triangular cell whose
third axis is periodic, with (hpx, hpy) = h₁ + h₂ the doubled in-plane
position and kp = k₁ + k₂ the doubled layer index:

    σ(hx, hy, k) = (mod(hpx − hx, 2 nx), mod(hpy − hy, 2 ny), mod(kp − k, nz)).

The in-plane part is the single-layer map [`_tri_reflect_site`](@ref);
the layer index is reflected alongside. The site set of the stacked cell,
{P(h) + k a₃}, is the Bravais lattice generated by the two primitive
translations of the triangular layer, (a, 0) and (a/2, √3 a/2), and the
third lattice vector a₃ = (Δ, h) (the two-site basis being the
centred-rectangular pair (0, 0), (1/2, √3/2)), so it is centrosymmetric
about every site and every two-site midpoint whatever the in-plane
stacking offset Δ (aligned layers, Δ = 0, or offset layers with Δ inside
the triangle). The map contains no Δ, so one method covers both
stackings. The layer wrap is exact because the periodic third supercell
vector is nz a₃ for every nz and every Δ: the site (h, k + nz) is the
periodic image of the site (h, k). Any pivot sends layer k to kp − k
(mod nz), an involution on layers; a site pivot (kp = 2 k₁, hpx even)
preserves the basis index and fixes its own layer (and, at even nz, the
layer k₁ + nz/2 as well), sending layer k to 2 k₁ − k (mod nz); at
nz = 2 an odd kp exchanges the two layers, and at nz = 3 any pivot pairs
the two layers whose indices sum to kp (mod 3) and fixes the third. σ is
an involution on site indices. At nz = 1 the layer term is identically
zero and σ reduces to the single-layer map.
"""
@inline function _tri_reflect_site_3d(site::Int, hpx::Int, hpy::Int, kp::Int,
                                      nx::Int, ny::Int, nz::Int)
    hx, hy, k = _tri_site_to_halfgrid(site, nx, ny)
    rx = mod(hpx - hx, 2 * nx)
    ry = mod(hpy - hy, 2 * ny)
    rk = mod(kp - k, nz)
    return _tri_halfgrid_to_site(rx, ry, rk, nx, ny)
end

"""
    geometric_cluster_swap!(lattice::MLattice{C,TriangularLattice}, p::Float64) where C

Geometric cluster move on the two-site-basis triangular lattice, single
layer or stacked, by point inversion through a random centrosymmetry
center.

The pivot is the midpoint of two random sites drawn with replacement,
which covers exactly the two families of inversion centers of the
triangular lattice (sites and half-lattice midpoints) and can never
produce an invalid pivot: plaquette centers, whose doubled position is
not a lattice vector, are unreachable by construction. The reflection
runs in integer half-grid coordinates ((hx, hy) = (2·ci + b, 2·cj + b) on
the (2·nx, 2·ny) torus), so no floating-point rounding or lookup table is
involved. Cluster growth and the pair swap are shared with the square
method via [`_build_geometric_cluster`](@ref).

On a stacked cell (`supercell_dimensions[3] > 1`, the layered-host
geometry of triangular layers stacked along a periodic third axis, either
aligned or with an in-plane stacking offset) the reflection is the
three-dimensional point inversion [`_tri_reflect_site_3d`](@ref): the
pivot's doubled layer index is the sum of the two pivot sites' layer
indices, and each site's layer is reflected alongside its in-plane
position. The stacked site set is a Bravais lattice for any stacking
offset, so the inversion is a symmetry of the whole site set and the
first-shell cluster growth may cross layers freely. Offset stacks have a
non-orthogonal supercell and must be built with `image_multiplicity=true`:
the default minimum-image neighbour path takes the rounded image rather
than the nearest one on a non-orthogonal supercell and can drop in-cutoff
bonds (pre-existing neighbour-list behaviour, not changed here). At
`supercell_dimensions[3] == 1` the layer term is identically zero, the
pivot and seed draws are the same draws in the same order, and the
integer arithmetic is the same, so single-layer trajectories are
unchanged.

Detailed balance holds exactly as in the square case: the reflection is
an involution on site indices, cluster growth is
configuration-independent, and the swap is symmetric, so the proposal
needs no Metropolis correction and nested sampling keeps its bare
energy-ceiling accept test (Heringa & Blöte, Phys. Rev. E 57, 4976
(1998)).

Requires the two-site basis, and on a stacked cell a periodic third axis
(`periodicity[3] == true`, since the layer reflection wraps on the c
period); violations throw an `ArgumentError`. A stacked cell with a
non-periodic third axis is not supported by this method.
"""
function geometric_cluster_swap!(lattice::MLattice{C,TriangularLattice}, p::Float64;
                                 record::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing) where C
    if length(lattice.basis) != 2
        throw(ArgumentError("geometric_cluster_swap! on a triangular lattice " *
            "requires the two-site centered-rectangular basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    nx, ny, nz = lattice.supercell_dimensions
    if nz != 1 && !lattice.periodicity[3]
        throw(ArgumentError("geometric_cluster_swap! on a stacked triangular " *
            "cell requires a periodic third axis (periodicity[3] == true), " *
            "since the layer reflection wraps on the c period; got $nz " *
            "layers with a non-periodic third axis"))
    end
    n_sites = num_sites(lattice)

    # Pivot: midpoint of two random sites (with replacement); keep its
    # doubled half-grid position and doubled layer index so all arithmetic
    # stays integer. At nz == 1 the layer indices are zero and the two
    # draws are those of the single-layer method.
    h1 = _tri_site_to_halfgrid(rand(1:n_sites), nx, ny)
    h2 = _tri_site_to_halfgrid(rand(1:n_sites), nx, ny)
    hpx = h1[1] + h2[1]
    hpy = h1[2] + h2[2]
    kp = h1[3] + h2[3]

    # Choose random seed site (1-indexed)
    seed = rand(1:n_sites)

    # Build cluster via BFS with fixed growth probability; the map reduces
    # to _tri_reflect_site at nz == 1
    reflect = s -> _tri_reflect_site_3d(s, hpx, hpy, kp, nx, ny, nz)
    cluster = _build_geometric_cluster(lattice, seed, reflect, p)

    # Swap occupation states for each applied (site, reflected_site) pair;
    # the recorded pairs revert the move when re-applied (involution)
    _apply_cluster_pairs!(lattice, cluster)
    if record !== nothing
        for pr in cluster
            pr[1] != pr[2] && push!(record, pr)
        end
    end

    return lattice
end

"""
    geometric_cluster_swap!(lattice::MLattice{C,GenericLattice}, p::Float64)

Guard method: geometric cluster moves are defined only for the square and triangular
geometries, whose reflection maps are integer grid involutions on the site indices; a
`GenericLattice` carries no such map, so this method always throws a descriptive
`ArgumentError` before drawing any randomness. Use swap-only decorrelation
(`clusters_freq = 0`) for generic geometries.
"""
function geometric_cluster_swap!(lattice::MLattice{C,GenericLattice}, p::Float64;
                                 record::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing) where C
    throw(ArgumentError("geometric_cluster_swap! supports square and triangular " *
        "lattices only, got a GenericLattice configuration; use swap-only " *
        "decorrelation (clusters_freq = 0) for generic geometries"))
end
