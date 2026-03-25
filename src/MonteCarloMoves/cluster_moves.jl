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

# Returns
- `lattice::MLattice{C,SquareLattice}`: The mutated lattice after the cluster swap.

# Notes
- Preserves per-component particle counts exactly.
- The move is self-inverse for a fixed cluster: applying the same cluster swap twice
  restores the original configuration.
- Currently restricted to `SquareLattice` with a single basis atom and
  `supercell_dimensions[3] == 1`.

# References
- Heringa & Blöte, Phys. Rev. E 57, 4976 (1998) — geometric cluster MC framework.
- Adaptation uses fixed growth probability (configuration-independent) for symmetric proposals.
"""
function geometric_cluster_swap!(lattice::MLattice{C,SquareLattice}, p::Float64) where C
    Lx = lattice.supercell_dimensions[1]
    Ly = lattice.supercell_dimensions[2]
    n_sites = num_sites(lattice)

    # Choose random pivot (grid coordinates, 0-indexed)
    pivot_gx = rand(0:Lx-1)
    pivot_gy = rand(0:Ly-1)

    # Choose random seed site (1-indexed)
    seed = rand(1:n_sites)

    # Build cluster via BFS with fixed growth probability
    cluster = _build_geometric_cluster(lattice, seed, pivot_gx, pivot_gy, Lx, Ly, p)

    # Swap occupation states for each (site, reflected_site) pair
    for (a, b) in cluster
        if a != b
            for comp in 1:C
                lattice.components[comp][a], lattice.components[comp][b] =
                    lattice.components[comp][b], lattice.components[comp][a]
            end
        end
    end

    return lattice
end

# --- Internal helpers (not exported) ---

"""
    _site_to_grid(site::Int, Lx::Int) -> Tuple{Int,Int}

Convert a 1-indexed site index to 0-indexed grid coordinates (gx, gy)
for a single-basis square lattice with Lx sites along x.
"""
@inline function _site_to_grid(site::Int, Lx::Int)
    s = site - 1
    gx = s % Lx
    gy = s ÷ Lx
    return (gx, gy)
end

"""
    _grid_to_site(gx::Int, gy::Int, Lx::Int) -> Int

Convert 0-indexed grid coordinates to a 1-indexed site index.
"""
@inline function _grid_to_site(gx::Int, gy::Int, Lx::Int)
    return gy * Lx + gx + 1
end

"""
    _reflect_site(site::Int, pivot_gx::Int, pivot_gy::Int, Lx::Int, Ly::Int) -> Int

Compute the point inversion of `site` through the pivot at grid coordinates
`(pivot_gx, pivot_gy)` with periodic wrapping on an Lx × Ly lattice.
"""
@inline function _reflect_site(site::Int, pivot_gx::Int, pivot_gy::Int, Lx::Int, Ly::Int)
    gx, gy = _site_to_grid(site, Lx)
    rx = mod(2 * pivot_gx - gx, Lx)
    ry = mod(2 * pivot_gy - gy, Ly)
    return _grid_to_site(rx, ry, Lx)
end

"""
    _build_geometric_cluster(lattice, seed, pivot_gx, pivot_gy, Lx, Ly, p) -> Vector{Tuple{Int,Int}}

Grow a geometric cluster from `seed` using BFS with fixed growth probability `p`.
Returns a vector of (site, reflected_site) pairs.
"""
function _build_geometric_cluster(lattice::MLattice{C,SquareLattice},
                                  seed::Int,
                                  pivot_gx::Int, pivot_gy::Int,
                                  Lx::Int, Ly::Int,
                                  p::Float64) where C
    visited = Set{Int}()
    stack = Int[seed]
    cluster = Tuple{Int,Int}[]

    while !isempty(stack)
        s = pop!(stack)

        s in visited && continue

        push!(visited, s)
        r = _reflect_site(s, pivot_gx, pivot_gy, Lx, Ly)
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
