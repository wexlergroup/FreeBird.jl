"""
    pbc_displacement(pos1, pos2, at::AbstractSystem)

Minimum-image displacement vector from `pos2` to `pos1` (`pos1 - pos2`,
min-imaged per axis), sharing `pbc_dist`'s per-axis periodicity switch and its
orthorhombic scope; its norm matches `pbc_dist(pos1, pos2, at)` on every axis
combination. Positions are assumed wrapped into the cell, the storage
convention of the walk kernels.

# Arguments
- `pos1`, `pos2`: The two positions (Unitful, Å components).
- `at::AbstractSystem`: The system carrying the cell and boundary conditions.

# Returns
- `SVector{3}`: The minimum-image displacement vector, in Å.
"""
function pbc_displacement(pos1, pos2, at::AbstractSystem)
    cellv = cell_vectors(at)
    pbc = periodicity(at)
    δx = pos1[1] - pos2[1]
    δy = pos1[2] - pos2[2]
    δz = pos1[3] - pos2[3]
    if pbc[1]
        L = cellv[1][1]
        δx -= L * round(δx / L)
    end
    if pbc[2]
        L = cellv[2][2]
        δy -= L * round(δy / L)
    end
    if pbc[3]
        L = cellv[3][3]
        δz -= L * round(δz / L)
    end
    return SVector(δx, δy, δz)
end

"""
    pair_force(r::typeof(1.0u"Å"), pot::SingleComponentPotential{Pairwise})
    pair_force(r::typeof(1.0u"Å"), lj::LJParameters)

The scalar radial force -du/dr of a pairwise potential at separation `r`, in
eV/Å. The generic method is a two-sided central finite difference of
`pair_energy` at a machine-scaled step, safe for any smooth pairwise potential;
`LJParameters` carries the analytic override
F(r) = 24 ε (2 (σ/r)¹² − (σ/r)⁶) / r for r ≤ r_c σ and exactly zero beyond:
plainly truncated, so the force is discontinuous at the cutoff (the stored
`shift` field shifts the energy only and never enters the derivative). In an
indicator-acceptance reflective walk the force enters only through the
reflection direction, where any deterministic position-dependent direction
preserves the stationary measure, so the finite-difference fallback degrades
mixing at worst, never sampled averages.
"""
function pair_force(r::typeof(1.0u"Å"), pot::SingleComponentPotential{Pairwise})
    h = (cbrt(eps(Float64)) * max(ustrip(u"Å", r), 1.0)) * u"Å"
    return (pair_energy(r - h, pot) - pair_energy(r + h, pot)) / (2 * h)
end

function pair_force(r::typeof(1.0u"Å"), lj::LJParameters)
    r > lj.cutoff * lj.sigma && return 0.0u"eV/Å"
    r6 = (lj.sigma / r)^6
    r12 = r6^2
    return 24 * lj.epsilon * (2 * r12 - r6) / r
end

"""
    interacting_gradient(at::AbstractSystem, pot::SingleComponentPotential{Pairwise},
                         list_num_par::Vector{Int}, frozen::Vector{Bool})

The per-particle gradient of the interacting energy, `∇ᵢU`, as a vector of
`SVector{3}` in eV/Å, over the same pair split as the matching
`interacting_energy` method: free-free and free-frozen pairs contribute,
frozen-frozen pairs are skipped, and frozen particles carry exactly zero
entries (they are never displaced; they still exert forces on free particles).
O(N²) serial accumulation through `pbc_dist`-consistent minimum-image
displacements; a collective move consumes the full gradient, so there is no
single-site variant to prefer. Empty and single-particle systems short-circuit
to zero entries.
"""
function interacting_gradient(at::AbstractSystem, pot::SingleComponentPotential{Pairwise},
                              list_num_par::Vector{Int}, frozen::Vector{Bool})
    n = length(at)
    zero_g = SVector(0.0u"eV/Å", 0.0u"eV/Å", 0.0u"eV/Å")
    grad = fill(zero_g, n)
    n <= 1 && return grad
    # Per-particle frozen status from the per-component flags
    bounds = cumsum(list_num_par)
    comp_of = [searchsortedfirst(bounds, i) for i in 1:n]
    part_frozen = [frozen[comp_of[i]] for i in 1:n]
    for i in 1:(n - 1)
        pos_i = position(at, i)
        for j in (i + 1):n
            (part_frozen[i] && part_frozen[j]) && continue
            disp = pbc_displacement(pos_i, position(at, j), at)
            r = sqrt(disp[1]^2 + disp[2]^2 + disp[3]^2)
            f = pair_force(r, pot)
            iszero(ustrip(f)) && continue
            # ∇ᵢU = (du/dr) r̂ᵢⱼ = -F(r) (posᵢ - posⱼ)/r
            g = (-f / r) * disp
            part_frozen[i] || (grad[i] += g)
            part_frozen[j] || (grad[j] -= g)
        end
    end
    return grad
end
