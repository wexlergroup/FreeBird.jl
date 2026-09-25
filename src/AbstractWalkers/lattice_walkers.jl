"""Whether independent fractional-coordinate wrapping is exact for this cell."""
@inline function _orthogonal_periodic_axes(
        cell::AbstractMatrix{Float64}, periodicity::Tuple{Bool,Bool,Bool})
    scale = maximum(norm(view(cell, :, k)) for k in 1:3)
    tol = 32 * eps(Float64) * max(scale^2, 1.0)
    for i in 1:3
        periodicity[i] || continue
        for j in 1:3
            i == j && continue
            norm(view(cell, :, j)) == 0.0 && continue
            abs(dot(view(cell, :, i), view(cell, :, j))) <= tol || return false
        end
    end
    return true
end

"""
    _minimum_image_displacement(cell, reciprocal, periodicity, dr)

Closest periodic image of Cartesian displacement `dr`. Orthogonal periodic
axes use direct fractional wrapping; skewed cells search every lattice
translation that can improve the wrapped candidate.
"""
@inline function _minimum_image_displacement(
        supercell_lattice_vectors::AbstractMatrix{Float64},
        reciprocal_lattice_vectors::AbstractMatrix{Float64},
        periodicity::Tuple{Bool,Bool,Bool}, dr)
    fractional_dr = reciprocal_lattice_vectors * dr
    center = ntuple(3) do k
        periodicity[k] ? -round(Int, fractional_dr[k]) : 0
    end

    n1, n2, n3 = center
    best_x = dr[1] + supercell_lattice_vectors[1, 1] * n1 +
                     supercell_lattice_vectors[1, 2] * n2 +
                     supercell_lattice_vectors[1, 3] * n3
    best_y = dr[2] + supercell_lattice_vectors[2, 1] * n1 +
                     supercell_lattice_vectors[2, 2] * n2 +
                     supercell_lattice_vectors[2, 3] * n3
    best_z = dr[3] + supercell_lattice_vectors[3, 1] * n1 +
                     supercell_lattice_vectors[3, 2] * n2 +
                     supercell_lattice_vectors[3, 3] * n3
    _orthogonal_periodic_axes(supercell_lattice_vectors, periodicity) &&
        return [best_x, best_y, best_z]
    best_distance2 = best_x^2 + best_y^2 + best_z^2
    best_distance = sqrt(best_distance2)

    # If another image is closer than `best`, its fractional component k has
    # magnitude at most ‖row_k(A⁻¹)‖ * ‖best‖. These bounds therefore enumerate
    # every lattice translation that can improve the result, including on
    # skewed cells where independently rounding fractional coordinates is not
    # a closest-vector algorithm.
    ranges = ntuple(3) do k
        if periodicity[k]
            bound = norm(view(reciprocal_lattice_vectors, k, :)) * best_distance
            slack = 32 * eps(Float64) * max(abs(fractional_dr[k]) + bound, 1.0)
            lo = min(center[k], ceil(Int, -fractional_dr[k] - bound - slack))
            hi = max(center[k], floor(Int, -fractional_dr[k] + bound + slack))
            lo:hi
        else
            0:0
        end
    end
    for n1 in ranges[1], n2 in ranges[2], n3 in ranges[3]
        x = dr[1] + supercell_lattice_vectors[1, 1] * n1 +
                    supercell_lattice_vectors[1, 2] * n2 +
                    supercell_lattice_vectors[1, 3] * n3
        y = dr[2] + supercell_lattice_vectors[2, 1] * n1 +
                    supercell_lattice_vectors[2, 2] * n2 +
                    supercell_lattice_vectors[2, 3] * n3
        z = dr[3] + supercell_lattice_vectors[3, 1] * n1 +
                    supercell_lattice_vectors[3, 2] * n2 +
                    supercell_lattice_vectors[3, 3] * n3
        distance2 = x^2 + y^2 + z^2
        if distance2 < best_distance2
            best_x, best_y, best_z = x, y, z
            best_distance2 = distance2
        end
    end
    return [best_x, best_y, best_z]
end

"""Minimum-image distance under the supplied cell and periodicity."""

@inline function _minimum_image_distance(supercell_lattice_vectors::AbstractMatrix{Float64},
                                         reciprocal_lattice_vectors::AbstractMatrix{Float64},
                                         periodicity::Tuple{Bool,Bool,Bool},
                                         pos_i, pos_j)
    dr = [pos_j[1] - pos_i[1], pos_j[2] - pos_i[2], pos_j[3] - pos_i[3]]
    return norm(_minimum_image_displacement(
        supercell_lattice_vectors, reciprocal_lattice_vectors, periodicity, dr))
end

"""
    compute_neighbors(supercell_lattice_vectors::Matrix{Float64},
                      positions::Matrix{Float64},
                      periodicity::Tuple{Bool, Bool, Bool},
                      cutoff_radii::Vector{Float64};
                      image_multiplicity::Bool=false)

Compute the neighbor shells of every site in a supercell. Each in-cutoff
distance is assigned to the first shell whose cutoff admits it (a nested
`<=` cutoff ladder), so `cutoff_radii` must be non-empty, finite, strictly
positive, and strictly increasing; anything else throws an
`ArgumentError`.

# Arguments
- `supercell_lattice_vectors::Matrix{Float64}`: The lattice vectors of the supercell, one cell vector per column.
- `positions::Matrix{Float64}`: The Cartesian positions of the sites, one site per row.
- `periodicity::Tuple{Bool, Bool, Bool}`: A Boolean tuple of length three indicating periodicity in each dimension (true for periodic, false for non-periodic).
- `cutoff_radii::Vector{Float64}`: The cutoff radii for the *index*-th nearest neighbors, strictly increasing.
- `image_multiplicity::Bool=false`: The neighbor-counting convention, see below.

# Returns
- `neighbors::Vector{Vector{Vector{Int}}}`: For each site, one vector of neighbor indices per shell.

# Conventions

With `image_multiplicity=false` (the default), each site pair is counted
once, in the shell of its minimum-image distance. On a cell whose periodic
circumference does not exceed twice a shell cutoff, pairs connected
through more than one periodic image are still counted once — the affected
shells sit below their bulk-tiled coordination — and one warning listing
every such shell and its collapsed-image count is emitted per call.
Detection is exact (the periodic images are enumerated), so faithful cells
never warn.

With `image_multiplicity=true`, a neighbor index is pushed once per
in-cutoff periodic image, including a site's own images (self-entries,
`j == i`, which arise when a periodic circumference is within the cutoff):
the cluster-expansion small-cell convention, under which the periodic
cell's energy equals the bulk energy per cell of the tiled configuration.
"""
function compute_neighbors(supercell_lattice_vectors::Matrix{Float64},
                           positions::Matrix{Float64},
                           periodicity::Tuple{Bool, Bool, Bool},
                           cutoff_radii::Vector{Float64};
                           image_multiplicity::Bool=false,
                           )

    if isempty(cutoff_radii) || any(r -> !(isfinite(r) && r > 0.0), cutoff_radii) ||
       !all(cutoff_radii[k] < cutoff_radii[k+1] for k in 1:length(cutoff_radii)-1)
        throw(ArgumentError(
            "cutoff_radii must be finite, positive, and strictly increasing " *
            "(shells are assigned by a nested cutoff ladder), got $cutoff_radii"))
    end

    neighbors = Vector{Vector{Vector{Int}}}(undef, size(positions, 1))
    num_atoms = size(positions, 1)

    # Compute reciprocal lattice vectors for minimum image convention
    a1 = supercell_lattice_vectors[:, 1]
    a2 = supercell_lattice_vectors[:, 2]
    a3 = supercell_lattice_vectors[:, 3]
    # A genuinely two-dimensional cell may carry a zero third lattice vector.
    # In that case the full 3x3 cell is singular even though its in-plane
    # lattice is valid. Build the reciprocal map from the 2x2 in-plane block;
    # the non-periodic Cartesian z displacement is restored after the
    # minimum-image round trip below.
    singular_2d = !periodicity[3] && all(iszero, a3)
    if singular_2d
        reciprocal_lattice_vectors = zeros(3, 3)
        reciprocal_lattice_vectors[1:2, 1:2] =
            inv([a1[1:2] a2[1:2]])
    else
        reciprocal_lattice_vectors = inv([a1 a2 a3])
    end

    layers_of_neighbors = length(cutoff_radii)
    r_max = last(cutoff_radii)
    orthogonal_periodic_axes =
        _orthogonal_periodic_axes(supercell_lattice_vectors, periodicity)

    # First shell whose cutoff admits the distance (the nested `<=` ladder);
    # 0 when the distance is beyond the outermost cutoff.
    function shell_of(distance::Float64)
        for k in 1:layers_of_neighbors
            if distance <= cutoff_radii[k]
                return k
            end
        end
        return 0
    end

    # If an image of displacement `dr` lies within r_max, its fractional
    # coordinate k lies within r_max * ||row_k(A^-1)|| of zero. The per-pair
    # ranges below center that bound on the *original* fractional displacement;
    # centering a fixed symmetric range on a minimum image is not sufficient
    # for strongly skewed cells.
    image_bounds = ntuple(k ->
        r_max * norm(view(reciprocal_lattice_vectors, k, :)), 3)

    # Discarded in-cutoff images per shell (default mode only)
    collapsed = zeros(Int, layers_of_neighbors)

    for i in 1:num_atoms
        nth_neighbors = Vector{Int}[]
        for _ in 1:layers_of_neighbors
            push!(nth_neighbors, Int[])
        end

        for j in 1:num_atoms
            dr = [positions[j, 1] - positions[i, 1],
                  positions[j, 2] - positions[i, 2],
                  positions[j, 3] - positions[i, 3]]
            original_fractional_dr = reciprocal_lattice_vectors * dr
            if orthogonal_periodic_axes
                fractional_dr = reciprocal_lattice_vectors * dr
                for k in 1:3
                    periodicity[k] &&
                        (fractional_dr[k] -= round(fractional_dr[k]))
                end
                dr_c = supercell_lattice_vectors * fractional_dr
                singular_2d && (dr_c[3] = dr[3])
            else
                dr_c = _minimum_image_displacement(
                    supercell_lattice_vectors, reciprocal_lattice_vectors,
                    periodicity, dr)
            end

            if image_multiplicity
                # Push j once per in-cutoff periodic image (self-images
                # included); only the single zero-displacement self term is
                # skipped — it is not a bond.
                ranges = ntuple(3) do k
                    if periodicity[k]
                        slack = 32 * eps(Float64) * max(
                            abs(original_fractional_dr[k]) + image_bounds[k], 1.0)
                        lo = ceil(Int,
                            -original_fractional_dr[k] - image_bounds[k] - slack)
                        hi = floor(Int,
                            -original_fractional_dr[k] + image_bounds[k] + slack)
                        lo:hi
                    else
                        0:0
                    end
                end
                for n1 in ranges[1], n2 in ranges[2], n3 in ranges[3]
                    (j == i && n1 == 0 && n2 == 0 && n3 == 0) && continue
                    dx = dr[1] + a1[1] * n1 + a2[1] * n2 + a3[1] * n3
                    dy = dr[2] + a1[2] * n1 + a2[2] * n2 + a3[2] * n3
                    dz = dr[3] + a1[3] * n1 + a2[3] * n2 + a3[3] * n3
                    d = sqrt(dx^2 + dy^2 + dz^2)
                    k = shell_of(d)
                    k != 0 && push!(nth_neighbors[k], j)
                end
            else
                # Minimum-image convention: for j != i push j once, into the
                # shell of the central (minimum-image) distance, exactly as
                # before; tally every other in-cutoff image — the collapsed
                # images of j != i pairs and all in-cutoff self-images — per
                # that image's own shell.
                if j != i
                    kept_shell = shell_of(norm(dr_c))
                    kept_shell != 0 && push!(nth_neighbors[kept_shell], j)
                end
                ranges = ntuple(3) do k
                    if periodicity[k]
                        slack = 32 * eps(Float64) * max(
                            abs(original_fractional_dr[k]) + image_bounds[k], 1.0)
                        lo = ceil(Int,
                            -original_fractional_dr[k] - image_bounds[k] - slack)
                        hi = floor(Int,
                            -original_fractional_dr[k] + image_bounds[k] + slack)
                        lo:hi
                    else
                        0:0
                    end
                end
                for n1 in ranges[1], n2 in ranges[2], n3 in ranges[3]
                    (j == i && n1 == 0 && n2 == 0 && n3 == 0) && continue
                    dx = dr[1] + a1[1] * n1 + a2[1] * n2 + a3[1] * n3
                    dy = dr[2] + a1[2] * n1 + a2[2] * n2 + a3[2] * n3
                    dz = dr[3] + a1[3] * n1 + a2[3] * n2 + a3[3] * n3
                    d = sqrt(dx^2 + dy^2 + dz^2)
                    k = shell_of(d)
                    k != 0 && (collapsed[k] += 1)
                end
                # One in-cutoff image of a distinct site is retained by the
                # minimum-image list rather than collapsed. Subtract it after
                # counting all images; this also handles equal-distance ties.
                j != i && kept_shell != 0 && (collapsed[kept_shell] -= 1)
            end
        end

        neighbors[i] = nth_neighbors
    end

    if any(>(0), collapsed)
        wrapped = ["$k (cutoff $(cutoff_radii[k])): $(collapsed[k]) collapsed image bond(s)"
                   for k in 1:layers_of_neighbors if collapsed[k] > 0]
        @warn "neighbor shell(s) [" * join(wrapped, ", ") * "] wrap the " *
              "periodic cell: some pairs (or a site and its own periodic " *
              "image) are connected through more than one image within the " *
              "cutoff but are counted once under the minimum-image " *
              "convention, so per-site coordination, and any energy " *
              "coupling these shells, is below the bulk-tiled value. " *
              "Enlarge the supercell so every periodic circumference " *
              "exceeds twice the shell cutoff, or construct the lattice " *
              "with `image_multiplicity=true` to count every image (the " *
              "cluster-expansion small-cell convention)."
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

_atomic_positions(slab) =
    [pyconvert(Vector{Float64}, atom.position) for atom in slab]

# ASE's named elemental surface builders.  The second entry is the FreeBird
# geometry of a translational adsorption-site orbit on that surface.  The
# stepped/rectangular faces deliberately use GenericLattice: calling them
# square merely because ASE happens to return an orthogonal cell would give
# them square-lattice order parameters that have no physical meaning.
const _ATOMIC_SURFACE_GEOMETRIES = Dict{Symbol,Symbol}(
    :fcc100     => :SquareLattice,
    :fcc110     => :GenericLattice,
    :fcc111     => :TriangularLattice,
    :fcc211     => :GenericLattice,
    :bcc100     => :SquareLattice,
    :bcc110     => :GenericLattice,
    :bcc111     => :TriangularLattice,
    :hcp0001    => :TriangularLattice,
    :hcp10m10   => :GenericLattice,
    :diamond100 => :SquareLattice,
    :diamond111 => :TriangularLattice,
)

"""Normalize and validate an ASE elemental surface-builder name."""
function _atomic_surface(surface::Union{Symbol,AbstractString})
    normalized = Symbol(lowercase(replace(string(surface), r"[()_\-]" => "")))
    normalized == :hcp1010 && (normalized = :hcp10m10)
    haskey(_ATOMIC_SURFACE_GEOMETRIES, normalized) || throw(ArgumentError(
        "unsupported ASE surface '$surface'; expected one of " *
        join(sort!(string.(collect(keys(_ATOMIC_SURFACE_GEOMETRIES)))), ", ")))
    return normalized
end

"""Construct one of ASE's named elemental slabs."""
function _build_atomic_surface(surface::Symbol, lattice_atom::String,
                               dimensions::Tuple{Int64,Int64,Int64},
                               lattice_constant::Float64,
                               lattice_constant_c::Union{Nothing,Float64}=nothing)
    builder = ase.build.__getattribute__(string(surface))
    if lattice_constant_c === nothing
        return builder(lattice_atom, dimensions; a=lattice_constant)
    end
    return builder(lattice_atom, dimensions;
                   a=lattice_constant, c=lattice_constant_c)
end

"""
Return the symmetry-equivalent offsets of an ASE named adsorption site.

ASE stores one representative position per named site.  Square-surface bridge
sites have two rotationally equivalent orientations and close-packed
triangular surfaces have three; all other names in the supported builders are
already distinct translational orbits (for example `longbridge` versus
`shortbridge`, and `fcc` versus `hcp`).
"""
function _surface_site_offsets(surface::Symbol, site::String,
                               representative::Tuple{Float64,Float64})
    if site == "bridge" && surface in (:fcc100, :bcc100)
        return [(0.5, 0.0), (0.0, 0.5)]
    elseif site == "bridge" && surface in (:fcc111, :hcp0001)
        return [(0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
    end
    return [representative]
end

"""
Build fcc(100) adsorption sites for periodic, finite, or mixed in-plane
boundary conditions.

ASE supplies one fractional representative for each named site family. Bridge
sites have two rotationally equivalent representatives. A representative with
a fractional offset along a finite direction is omitted from the final unit
cell in that direction, while a periodic direction retains the boundary-crossing
site. This gives `nx*ny` ontop sites,
`(nx - !px)*ny + nx*(ny - !py)` bridge sites, and
`(nx - !px)*(ny - !py)` hollow sites.
"""
function _fcc100_surface_sites(slab,
        dimensions::Tuple{Int64,Int64,Int64},
        periodicity::Tuple{Bool,Bool,Bool},
        type_of_sites::Vector{String})
    info = slab.info["adsorbate_info"]
    named = info["sites"]
    available = sort!(pyconvert(Vector{String}, pylist(named.keys())))
    invalid = filter(site -> site ∉ available, type_of_sites)
    isempty(invalid) || throw(ArgumentError(
        "adsorption site type(s) $(invalid) are not available on fcc100; " *
        "ASE provides $(available)"))

    unit_cell = pyconvert(Matrix{Float64}, info["cell"])
    cell = pyconvert(Matrix{Float64}, slab.get_cell())
    nx, ny, _ = dimensions
    sites = Tuple{Float64,Float64}[]
    for site in type_of_sites
        raw = pyconvert(Vector{Float64}, named[site])
        representative = (raw[1], raw[2])
        for offset in _surface_site_offsets(:fcc100, site, representative)
            crosses_x = !isapprox(mod(offset[1], 1.0), 0.0; atol=1e-12)
            crosses_y = !isapprox(mod(offset[2], 1.0), 0.0; atol=1e-12)
            ni = nx - (!periodicity[1] && crosses_x)
            nj = ny - (!periodicity[2] && crosses_y)
            (ni == 0 || nj == 0) && continue
            for j in 0:(nj - 1), i in 0:(ni - 1)
                u, v = i + offset[1], j + offset[2]
                push!(sites,
                      (u * unit_cell[1, 1] + v * unit_cell[2, 1],
                       u * unit_cell[1, 2] + v * unit_cell[2, 2]))
            end
        end
    end
    return _fold_surface_sites(sites, cell)
end

"""Fold and de-duplicate Cartesian in-plane sites in an ASE slab cell."""
function _fold_surface_sites(sites::Vector{Tuple{Float64,Float64}},
                             cell::Matrix{Float64}; tol::Float64=1e-8)
    inplane = [cell[1, 1] cell[2, 1]; cell[1, 2] cell[2, 2]]
    reciprocal = inv(inplane)
    fractional = Tuple{Float64,Float64}[]
    for (x, y) in sites
        f = mod.(reciprocal * [x, y], 1.0)
        any(g -> norm(f .- collect(g) .- round.(f .- collect(g))) <= tol,
            fractional) || push!(fractional, (f[1], f[2]))
    end
    return map(fractional) do f
        xy = inplane * collect(f)
        (xy[1], xy[2])
    end
end

"""
Build adsorption sites from the metadata supplied by ASE's surface builder.

Each selected name contributes its complete translational orbit.  Rotationally
equivalent `bridge` orientations are expanded explicitly.  `fcc211` is the one
ASE elemental builder without named sites, so FreeBird exposes its top-layer
atoms as `ontop` sites only.
"""
function _ase_surface_sites(slab, surface::Symbol,
                            dimensions::Tuple{Int64,Int64,Int64},
                            type_of_sites::Vector{String})
    positions = _atomic_positions(slab)
    cell = pyconvert(Matrix{Float64}, slab.get_cell())

    if surface == :fcc211
        type_of_sites == ["ontop"] || throw(ArgumentError(
            "ASE's fcc211 builder provides no named adsorption sites; " *
            "AtomicLattice supports type_of_sites=[\"ontop\"] for fcc211"))
        z_top = maximum(p[3] for p in positions)
        sites = [(p[1], p[2]) for p in positions if isapprox(p[3], z_top; atol=1e-8)]
        return _fold_surface_sites(sites, cell)
    end

    info = slab.info["adsorbate_info"]
    named = info["sites"]
    available = sort!(pyconvert(Vector{String}, pylist(named.keys())))
    invalid = filter(site -> site ∉ available, type_of_sites)
    isempty(invalid) || throw(ArgumentError(
        "adsorption site type(s) $(invalid) are not available on $surface; " *
        "ASE provides $(available)"))

    unit_cell = pyconvert(Matrix{Float64}, info["cell"])
    nx, ny, _ = dimensions
    sites = Tuple{Float64,Float64}[]
    for site in type_of_sites
        raw = pyconvert(Vector{Float64}, named[site])
        representative = (raw[1], raw[2])
        for offset in _surface_site_offsets(surface, site, representative)
            for j in 0:(ny - 1), i in 0:(nx - 1)
                u, v = i + offset[1], j + offset[2]
                push!(sites,
                      (u * unit_cell[1, 1] + v * unit_cell[2, 1],
                       u * unit_cell[1, 2] + v * unit_cell[2, 2]))
            end
        end
    end
    isempty(sites) && throw(ArgumentError(
        "type_of_sites must select at least one adsorption-site family"))
    return _fold_surface_sites(sites, cell)
end

"""Return nested cutoffs for the first `n_shells` periodic site distances."""
function _site_shell_cutoffs(positions::Matrix{Float64}, cell::Matrix{Float64},
                             periodicity::Tuple{Bool,Bool,Bool}, n_shells::Int;
                             image_multiplicity::Bool=false)
    n_shells == 0 && return Float64[]
    supercell = permutedims(cell) # ASE stores cell vectors as rows
    # ASE surface builders leave the non-periodic third cell vector at zero
    # unless the caller adds vacuum. Match `compute_neighbors`' 2-D reciprocal
    # construction so a perfectly valid slab is not rejected as singular.
    singular_2d = !periodicity[3] && all(iszero, view(supercell, :, 3))
    if singular_2d
        reciprocal = zeros(3, 3)
        reciprocal[1:2, 1:2] = inv(supercell[1:2, 1:2])
    else
        reciprocal = inv(supercell)
    end
    distances = Float64[]
    n_sites = size(positions, 1)
    if image_multiplicity
        radius = n_shells + 1
        ranges = ntuple(k -> periodicity[k] ? (-radius:radius) : (0:0), 3)
        for i in 1:n_sites, j in 1:n_sites,
            n1 in ranges[1], n2 in ranges[2], n3 in ranges[3]
            i == j && n1 == 0 && n2 == 0 && n3 == 0 && continue
            dx = positions[j, 1] - positions[i, 1] +
                 supercell[1, 1] * n1 + supercell[1, 2] * n2 + supercell[1, 3] * n3
            dy = positions[j, 2] - positions[i, 2] +
                 supercell[2, 1] * n1 + supercell[2, 2] * n2 + supercell[2, 3] * n3
            dz = positions[j, 3] - positions[i, 3] +
                 supercell[3, 1] * n1 + supercell[3, 2] * n2 + supercell[3, 3] * n3
            d = sqrt(dx^2 + dy^2 + dz^2)
            d > 1e-10 && push!(distances, d)
        end
    else
        n_sites >= 2 || throw(ArgumentError(
            "at least two adsorption sites are required to build neighbor shells " *
            "with image_multiplicity=false; use num_nearest_neighbors=0 for a " *
            "noninteracting/MLIP lattice or image_multiplicity=true to include " *
            "periodic self images"))
        for i in 1:(n_sites - 1), j in (i + 1):n_sites
            d = _minimum_image_distance(
                supercell, reciprocal, periodicity,
                view(positions, i, :), view(positions, j, :))
            d > 1e-10 && push!(distances, d)
        end
    end
    sort!(distances)
    shells = Float64[]
    for d in distances
        (isempty(shells) || !isapprox(d, last(shells); rtol=1e-8, atol=1e-10)) &&
            push!(shells, d)
    end
    length(shells) >= n_shells || throw(ArgumentError(
        "requested $n_shells adsorption-site neighbor shells, but this lattice " *
        "contains only $(length(shells)) distinct nonzero periodic distances"))
    return shells[1:n_shells] .* (1 + 1e-8)
end

"""
Build one site-index point-inversion map used by a geometric cluster move.

An empty row means the selected adsorption-site set is not closed under
periodic point inversion about `pivot`. A nonempty row is checked as both a
bijection and an involution. Building only the requested pivot keeps storage
linear in the number of sites.
"""
function _build_atomic_reflection_row(
        sites::Vector{Tuple{Float64,Float64}}, cell::Matrix{Float64},
        periodicity::Tuple{Bool,Bool,Bool}, tol::Float64, pivot::Int)
    periodicity[1] && periodicity[2] || return Int[]
    inplane = [cell[1, 1] cell[2, 1]; cell[1, 2] cell[2, 2]]
    inv_inplane = inv(inplane)
    fractional = [mod.(inv_inplane * [x, y], 1.0) for (x, y) in sites]
    n = length(sites)
    fractional_tol = tol / minimum(svdvals(inplane))
    n_bins = max(1, floor(Int, 1 / max(fractional_tol, eps(Float64))))
    bucket(f) = (mod(floor(Int, mod(f[1], 1.0) * n_bins), n_bins),
                 mod(floor(Int, mod(f[2], 1.0) * n_bins), n_bins))
    buckets = Dict{Tuple{Int,Int},Vector{Int}}()
    for (site, f) in enumerate(fractional)
        push!(get!(buckets, bucket(f), Int[]), site)
    end

    checkbounds(sites, pivot)
    row = zeros(Int, n)
    for site in 1:n
        reflected = mod.(2 .* fractional[pivot] .- fractional[site], 1.0)
        best = 0
        best_distance = Inf
        base = bucket(reflected)
        for di in -1:1, dj in -1:1
            key = (mod(base[1] + di, n_bins), mod(base[2] + dj, n_bins))
            for candidate in get(buckets, key, Int[])
                delta = fractional[candidate] .- reflected
                delta .-= round.(delta)
                distance = norm(inplane * delta)
                if distance < best_distance
                    best = candidate
                    best_distance = distance
                end
            end
        end
        best_distance <= tol || return Int[]
        row[site] = best
    end
    length(unique(row)) == n || return Int[]
    all(row[row[site]] == site for site in 1:n) || return Int[]
    return row
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
        cutoff_radii::Vector{Float64},
        components::Vector{Vector{Bool}},
        adsorptions::Vector{Bool};
        image_multiplicity::Bool=false,
    ) where {C,G}

Creates an `MLattice` instance with the specified parameters. The constructor performs the following steps:

1. Validates that the number of components matches the expected value `C`.
2. Computes the positions of the lattice points using `lattice_positions`.
3. Computes the supercell lattice vectors.
4. Computes the neighbors of each lattice point using `compute_neighbors`, passing the `image_multiplicity` keyword through.

Throws an `ArgumentError` if the number of components does not match `C`.

# Outer Constructors

    MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
                               interlayer_spacing::Union{Nothing,Float64}=nothing,
                               basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
                               supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
                               periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                               cutoff_radii::Vector{Float64}=[1.1, 1.5],
                               components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                               adsorptions::Union{Vector{Int},Symbol}=:full,
                               image_multiplicity::Bool=false)

    MLattice{C,TriangularLattice}(; lattice_constant::Float64=1.0,
                                  interlayer_spacing::Union{Nothing,Float64}=nothing,
                                  basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)],
                                  supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 2, 1),
                                  periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                  cutoff_radii::Vector{Float64}=[1.1, 1.8],
                                  components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                  adsorptions::Union{Vector{Int},Symbol}=:full,
                                  image_multiplicity::Bool=false,
                                  layer_offset::Union{Nothing,NTuple{2,Float64}}=nothing,
                                  stacking::Symbol=:aligned)

Constructs a square/triangular lattice with the specified parameters. The `components` and `adsorptions` arguments can be a vector of integers specifying
the indices of the occupied sites, or a symbol. If `components` is `:equal`, the lattice is divided into `C` equal components when possible, or
nearest to equal components otherwise. If `adsorptions` is `:full`, all sites are classified as adsorption sites.

The `image_multiplicity` keyword selects the neighbor-counting convention of
[`compute_neighbors`](@ref): under the default minimum-image convention each
site pair is counted once and a warning reports any shells that wrap the
periodic cell; with `image_multiplicity=true` every in-cutoff periodic image
is counted, self-image entries included (the cluster-expansion small-cell
convention, under which the cell's energy equals the bulk energy per cell of
the tiled configuration). The duplicated entries compose correctly
downstream: the energy kernels add `coupling/2` per ordered-pair entry, so a
doubled bond reaches the full bond energy and a self-image entry passes the
occupation test exactly when the site is occupied; geometric-cluster growth
reads only the first shell, and growth across a duplicated bond remains a
configuration-independent symmetric proposal.

The `interlayer_spacing` keyword sets the out-of-plane (third-axis) lattice
spacing. The default `nothing` means isotropic spacing: the third lattice
vector is `[0, 0, lattice_constant]`. Note this is a behavior change on one
previously broken path: 3D cells built with `lattice_constant ≠ 1` used to
mix scales (the out-of-plane spacing was fixed at 1.0, so a "nearest-neighbor"
cutoff could select only interlayer bonds without warning); such cells now
default to isotropic spacing. Pass `interlayer_spacing = 1.0` to reproduce
the old geometry. An explicit value `c` gives a tetragonal cell whose
square-lattice distance ladder is a, c, √2·a, √(a² + c²), 2a, …, so a
suitable cutoff ladder separates in-plane from interlayer nearest neighbors
into distinct shells: with a = 1.0, c = 1.25 and `cutoff_radii = [1.1, 1.35]`,
shell 1 is in-plane only and shell 2 interlayer only (the next distance
√2 ≈ 1.414 is excluded), so a `GenericLatticeHamiltonian` with couplings
`[J∥, J⊥]` expresses direction-resolved nearest-neighbor interactions with no
other library changes. Avoid degenerate spacings that alias the two bond
classes into one shell (c = √2·a puts interlayer bonds at the in-plane
second-neighbor distance, c = 2a at the third); the strictly-increasing
cutoff-ladder validation and the empty-shell warning of
[`compute_neighbors`](@ref) are the runtime backstops. The keyword composes
with `image_multiplicity`; the third axis is non-periodic by default, so
slabs gain no z-wrap images. Must be finite and strictly positive when given;
violations throw an `ArgumentError`.

The triangular constructor's `layer_offset` and `stacking` keywords set the
in-plane displacement between successive layers, for layered hosts whose
layers do not sit vertically above one another. The default
(`layer_offset = nothing`, `stacking = :aligned`) keeps the vertical third
lattice vector `[0, 0, c]`, so every existing call constructs the same
object. An explicit `layer_offset = (Δx, Δy)`, an absolute in-plane
displacement in the same length units as `lattice_constant` and
`interlayer_spacing`, makes the third lattice vector `[Δx, Δy, c]`, so
layer `k` sits at `(k − 1)·(Δx, Δy)` in-plane; `stacking = :aligned` with
an explicit offset is accepted. `stacking = :abc` selects the
triangle-centre offset `(a/2, √3·a/6)`: each layer sits over the triangle
centres of the layer below and the site geometry repeats after three
layers (`3·Δ = t₁ + t₂` is an in-plane lattice vector). It cannot be
combined with an explicit `layer_offset`, and any other `stacking` value
throws an `ArgumentError`. Both offset components must be finite. No
layer-count condition applies with a periodic third axis: the supercell's
third vector is `supercell_dimensions[3]` times the skewed `a₃`, itself a
lattice vector of the stacked crystal, so every layer count gives a valid
periodic crystal with the same ABC site geometry, and
`supercell_dimensions[3]` is the period of the occupation pattern along
the stacking direction: `1` constrains every layer to the same pattern in
its own frame (a cell primitive along the stacking direction), while `3`
or `6` admits the stacking sequences of an in-plane order (see
[`stacking_bragg_amplitude`](@ref)). On an offset cell build with
`image_multiplicity = true`: on one- and two-layer periodic cells the
default once-per-pair minimum-image convention drops the interlayer shell
entirely or by half, because the partners across the stacking direction
are further images of in-plane neighbours (a warning fires), and on
taller skewed cells it can still drop pairs from a long-cutoff shell,
whereas the all-images convention enumerates every in-cutoff image and
reproduces the bulk coordination of the stacked crystal. The skewed cell
needs no other change: the layer
helpers read the contiguous layer blocks of `lattice_positions`.

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
        adsorptions::Vector{Bool};
        image_multiplicity::Bool=false,
    ) where {C,G}

        num_components = length(components)

        if num_components != C
            throw(ArgumentError("For a $C-component system, got $num_components components!"))
        end

        positions = lattice_positions(lattice_vectors, basis, supercell_dimensions)

        supercell_lattice_vectors = lattice_vectors * Diagonal([supercell_dimensions[1], supercell_dimensions[2], supercell_dimensions[3]])
        neighbors = compute_neighbors(supercell_lattice_vectors, positions, periodicity, cutoff_radii;
                                      image_multiplicity=image_multiplicity)

        
        return new{C,G}(lattice_vectors, positions, basis, supercell_dimensions, periodicity, cutoff_radii, components, neighbors, adsorptions)
    end

    # Shared-geometry constructor: alias the eight run-invariant fields of
    # `source` (they are written only during construction; the Monte Carlo
    # kernels mutate occupancies exclusively) and install the fresh
    # `components`. Skips lattice_positions and the all-pairs
    # compute_neighbors entirely; see `replicate_walkers` for the exported
    # entry point.
    function MLattice{C,G}(::Val{:share_geometry},
                           source::MLattice{C,G},
                           components::Vector{Vector{Bool}}) where {C,G}
        if length(components) != C
            throw(ArgumentError("For a $C-component system, got $(length(components)) components!"))
        end
        return new{C,G}(source.lattice_vectors, source.positions,
                        source.basis, source.supercell_dimensions,
                        source.periodicity, source.cutoff_radii,
                        components, source.neighbors, source.adsorptions)
    end
end


"""
    mutable struct AtomicLattice{C,G} <: AbstractLattice

A mutable struct representing an atomic adsorption lattice using ASE
(Atomic Simulation Environment).

Multiple adsorbate species use the same mutually-exclusive component masks as
`MLattice`. The `surface` keyword selects one of ASE's named elemental surface
builders and determines whether `G` is `SquareLattice`, `TriangularLattice`, or
`GenericLattice`.

# Fields
- `lattice_atom::String`: The chemical symbol of the lattice substrate atom.
- `surface::Symbol`: ASE surface builder used for the substrate (for example
  `:fcc100`, `:fcc111`, or `:bcc110`).
- `adsorbate_atoms::Vector{String}`: The chemical symbols of the adsorbate species.
- `all_sites::Vector{Tuple{Float64, Float64}}`: Coordinates of every adsorption
  site, grouped by the requested site-family order and then by surface-cell
  order.
- `geometry_fingerprint::UInt64`: Constructor-time signature used to reject
  unsupported direct mutation of geometry-defining fields.
- `components::Vector{Vector{Bool}}`: **The ground truth.** One mask per
  adsorbate species, indexed over `all_sites`, with at most one species per site.
- `reflection_row::Vector{Int}`: Cached periodic point-inversion map for one
  cluster-move pivot.
- `reflection_pivot::Int`: Pivot represented by `reflection_row`.
- `reflection_row_built::Bool`: Whether the cached row has been checked.
- `adsorbate_height::Float64`: Height at which adsorbates are placed. Stored rather than assumed, so the constructor and `sync_ase_lattice!` cannot disagree about it.
- `ase_dirty::Bool`: Whether `ase_lattice` is stale with respect to `components`.
- `synced_components::Vector{Vector{Bool}}`: Occupancy snapshot represented by
  `ase_lattice`; it detects direct mutations of the public component masks.
- `synced_adsorbate_atoms::Vector{String}` and
  `synced_adsorbate_height::Float64`: adsorbate metadata represented by the
  cached ASE structure.
- `ase_cache_fingerprint::UInt64`: Signature of the synchronized ASE frame,
  used to reject direct mutation of the derived Python cache.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The dimensions of the supercell.
- `lattice_constant::Float64`: The lattice constant of the unit cell.
- `lattice_constant_c::Union{Nothing,Float64}`: Optional hcp c-axis lattice
  constant. `nothing` uses ASE's reference c/a ratio.
- `periodicity::Tuple{Bool, Bool, Bool}`: The periodic boundary conditions in each dimension.
- `lattice_positions::Matrix{Float64}`: The adsorption-site positions, with z=0.
- `num_nearest_neighbors::Int64`: The number of adsorption-site neighbor shells.
- `image_multiplicity::Bool`: Whether neighbor lists retain every in-cutoff
  periodic image, including self images, rather than one minimum-image entry
  per site pair.
- `cutoff_radii::Vector{Float64}`: Nested cutoffs for those shells.
- `neighbors::Vector{Vector{Vector{Int}}}`: Neighbor lists indexed over `all_sites`.
- `type_of_sites::Vector{String}`: The types of adsorption sites (e.g., "ontop", "bridge", "hollow").
- `ase_lattice::Py`: The ASE atoms object. A **derived cache** of `components`, not a second source of truth — see `sync_ase_lattice!`.

# Constructor
    AtomicLattice{C,G}(;
        lattice_atom::String,
        surface::Union{Symbol,AbstractString}=:fcc100,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        lattice_constant_c::Union{Nothing,Float64} = nothing,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String},
        coverage::Float64 = 0.5,
        components::Union{Nothing,AbstractVector{<:Integer},
                          AbstractVector{<:AbstractVector{Bool}}}=nothing,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String},
        adsorbate_height::Float64 = 1.0,
        image_multiplicity::Bool = false
    ) where {C,G}

Creates an `AtomicLattice` instance with the specified parameters. The constructor performs the following steps:
1. Validates the selected ASE surface, its matching FreeBird geometry, and its
   scalar and collection arguments.
2. Constructs the selected elemental slab using ASE with the specified lattice atom and dimensions.
3. Sets the periodic boundary conditions on the slab.
4. Adds adsorbates to the surface at the specified sites with the given coverage.
5. Computes the lattice positions and neighbor lists.

Throws an `ArgumentError` for an unsupported component/geometry type or invalid
constructor argument.

# Arguments
- `lattice_atom::String`: Chemical symbol for the substrate (e.g., "Pt", "Cu").
- `surface`: ASE elemental surface builder. Supported values are `fcc100`,
  `fcc110`, `fcc111`, `fcc211`, `bcc100`, `bcc110`, `bcc111`, `hcp0001`,
  `hcp10m10`, `diamond100`, and `diamond111`. Parentheses, underscores, and
  hyphens in string values are ignored.
- `supercell_dimensions::Tuple{Int64, Int64, Int64}`: The first two entries
  repeat the surface in plane; the third is the number of substrate layers
  passed to the ASE surface builder.
- `lattice_constant::Float64`: Lattice constant in Ångströms.
- `lattice_constant_c`: Optional hcp c-axis lattice constant in Ångströms.
  It is valid only for `hcp0001` and `hcp10m10`; `nothing` uses ASE's
  reference c/a ratio.
- `periodicity::Tuple{Bool, Bool, Bool}`: Periodic boundary conditions for each dimension.
- `adsorbate_atoms::Vector{String}`: Chemical symbols of adsorbates.
- `coverage::Float64`: Requested total fractional surface coverage (default:
  `0.5`). The constructor rounds this to the nearest realizable number of
  occupied sites. For multiple species, those sites are distributed as evenly
  as possible among the `C` species.
- `components`: Optional per-species particle counts, or explicit Boolean
  masks (including `BitVector`s) over the adsorption sites. When supplied,
  this overrides `coverage`.
- `num_nearest_neighbors::Int64`: Number of neighbor shells to construct. Use
  zero when the selected energy model and moves do not use neighbor lists.
- `type_of_sites::Vector{String}`: Adsorption-site families to include.
- `adsorbate_height::Float64`: Height in Ångströms used when adding every
  adsorbate to the ASE surface. This single construction height is shared by
  all adsorbate species; it is not relaxed automatically.
- `image_multiplicity::Bool`: Whether neighbor shells retain each in-cutoff
  periodic image separately. The default keeps one minimum-image entry per
  site pair.

# Returns
- `AtomicLattice{C,G}`: An atomic lattice object with `C` adsorbate species and geometry type `G`.
"""
mutable struct AtomicLattice{C,G} <: AbstractLattice
    lattice_atom::String
    surface::Symbol
    adsorbate_atoms::Vector{String}
    supercell_dimensions::Tuple{Int64, Int64, Int64}
    lattice_constant::Float64
    lattice_constant_c::Union{Nothing,Float64}
    periodicity::Tuple{Bool, Bool, Bool}
    lattice_positions::Matrix{Float64}
    num_nearest_neighbors::Int64
    image_multiplicity::Bool
    cutoff_radii::Vector{Float64}
    neighbors::Vector{Vector{Vector{Int}}}
    type_of_sites::Vector{String}
    all_sites::Vector{Tuple{Float64, Float64}}
    geometry_fingerprint::UInt64
    # ── ground truth and site-index geometry ───────────────────────────────
    components::Vector{Vector{Bool}}
    reflection_row::Vector{Int}
    reflection_pivot::Int
    reflection_row_built::Bool
    adsorbate_height::Float64
    # ── derived cache of the above; see sync_ase_lattice! ───────────────────
    ase_lattice::Py
    ase_dirty::Bool
    synced_components::Vector{Vector{Bool}}
    synced_adsorbate_atoms::Vector{String}
    synced_adsorbate_height::Float64
    ase_cache_fingerprint::UInt64

    function AtomicLattice{C,G}(;
        lattice_atom::String,
        surface::Union{Symbol,AbstractString}=:fcc100,
        supercell_dimensions::Tuple{Int64, Int64, Int64},
        lattice_constant::Float64,
        lattice_constant_c::Union{Nothing,Float64} = nothing,
        periodicity::Tuple{Bool, Bool, Bool},
        adsorbate_atoms::Vector{String},
        coverage::Float64 = 0.5,
        components::Union{Nothing,AbstractVector{<:Integer},
                          AbstractVector{<:AbstractVector{Bool}}}=nothing,
        num_nearest_neighbors::Int64,
        type_of_sites::Vector{String},
        adsorbate_height::Float64 = 1.0,
        image_multiplicity::Bool = false
    ) where {C,G}

        C > 0 || throw(ArgumentError("AtomicLattice requires at least one component"))
        surface = _atomic_surface(surface)
        expected_geometry_name = _ATOMIC_SURFACE_GEOMETRIES[surface]
        expected_geometry = getfield(@__MODULE__, expected_geometry_name)
        G === expected_geometry || throw(ArgumentError(
            "ASE surface $surface uses $expected_geometry_name in AtomicLattice, " *
            "but the requested geometry is $(nameof(G))"))
        !periodicity[3] || throw(ArgumentError(
            "AtomicLattice requires a non-periodic z direction because ASE " *
            "surface builders provide no out-of-plane cell vector; got " *
            "periodicity=$periodicity"))
        if surface != :fcc100 && !(periodicity[1] && periodicity[2])
            throw(ArgumentError(
                "AtomicLattice surface $surface requires periodic " *
                "x and y boundaries so ASE named adsorption sites form a " *
                "complete lattice; got periodicity=$periodicity"))
        end
        isfinite(lattice_constant) && lattice_constant > 0 || throw(ArgumentError(
            "lattice_constant must be finite and positive, got $lattice_constant"))
        if lattice_constant_c !== nothing
            surface in (:hcp0001, :hcp10m10) || throw(ArgumentError(
                "lattice_constant_c is supported only for hcp surfaces, got $surface"))
            isfinite(lattice_constant_c) && lattice_constant_c > 0 ||
                throw(ArgumentError("lattice_constant_c must be finite and positive, " *
                                    "got $lattice_constant_c"))
        end
        all(>(0), supercell_dimensions) || throw(ArgumentError(
            "supercell_dimensions must be positive, got $supercell_dimensions"))
        if surface == :fcc211 && supercell_dimensions[1] % 3 != 0
            throw(ArgumentError(
                "ASE fcc211 requires supercell_dimensions[1] divisible by 3; " *
                "got $(supercell_dimensions[1])"))
        elseif surface == :hcp10m10 && isodd(supercell_dimensions[2])
            throw(ArgumentError(
                "ASE hcp10m10 requires an even supercell_dimensions[2]; " *
                "got $(supercell_dimensions[2])"))
        end
        isfinite(coverage) && 0.0 <= coverage <= 1.0 || throw(ArgumentError(
            "coverage must be finite and between 0 and 1, got $coverage"))
        isfinite(adsorbate_height) || throw(ArgumentError(
            "adsorbate_height must be finite, got $adsorbate_height"))
        num_nearest_neighbors >= 0 || throw(ArgumentError(
            "num_nearest_neighbors must be nonnegative, got $num_nearest_neighbors"))

        num_adsorbates = length(adsorbate_atoms)

        if num_adsorbates != C
            throw(ArgumentError("For a $C-adsorbate system, got $num_adsorbates adsorbates"))
        end
        all(!isempty, adsorbate_atoms) || throw(ArgumentError(
            "adsorbate atom symbols must be non-empty"))
        isempty(type_of_sites) && throw(ArgumentError(
            "type_of_sites must select at least one adsorption-site family"))
        length(unique(type_of_sites)) == length(type_of_sites) || throw(ArgumentError(
            "type_of_sites may not contain duplicates, got $type_of_sites"))
        adsorbate_atoms = copy(adsorbate_atoms)
        type_of_sites = copy(type_of_sites)

        slab = _build_atomic_surface(
            surface, lattice_atom, supercell_dimensions, lattice_constant,
            lattice_constant_c)
        # Most ASE surface builders tag substrate layers with positive integers,
        # but fcc211 leaves every substrate atom at tag 0.  FreeBird reserves
        # tag 0 for atoms added by `add_adsorbate`, so normalize any untagged
        # substrate atoms before the ASE frame becomes a mutable cache.
        substrate_tags = pyconvert(Vector{Int}, slab.get_tags())
        any(==(0), substrate_tags) &&
            slab.set_tags([tag == 0 ? 1 : tag for tag in substrate_tags])
        slab.set_pbc(periodicity)

        if surface == :fcc100
            all_sites = _fcc100_surface_sites(
                slab, supercell_dimensions, periodicity, type_of_sites)
        else
            all_sites = _ase_surface_sites(
                slab, surface, supercell_dimensions, type_of_sites)
        end
        isempty(all_sites) && throw(ArgumentError(
            "no adsorption sites produced for type_of_sites=$type_of_sites " *
            "with dimensions=$supercell_dimensions and periodicity=$periodicity"))
        ase_lattice = slab
        occupied = falses(length(all_sites))
        if components === nothing
            n_occupied = round(Int, coverage * length(all_sites))
            if n_occupied == length(all_sites)
                occupied .= true
            elseif n_occupied > 0
                occupied[randperm(length(all_sites))[1:n_occupied]] .= true
            end
        end

        cell = pyconvert(Matrix{Float64}, slab.get_cell())
        lattice_positions = hcat(first.(all_sites), last.(all_sites),
                                 zeros(length(all_sites)))
        cutoff_radii = _site_shell_cutoffs(
            lattice_positions, cell, periodicity, num_nearest_neighbors;
            image_multiplicity)
        neighbors = if isempty(cutoff_radii)
            [Vector{Int}[] for _ in 1:length(all_sites)]
        else
            compute_neighbors(permutedims(cell), lattice_positions,
                              periodicity, cutoff_radii; image_multiplicity)
        end
        component_masks = [falses(length(all_sites)) for _ in 1:C]
        if components === nothing
            selected_sites = C == 1 ? findall(occupied) : shuffle(findall(occupied))
            for (k, site) in enumerate(selected_sites)
                component_masks[mod1(k, C)][site] = true
            end
        elseif components isa AbstractVector{<:Integer}
            length(components) == C || throw(ArgumentError(
                "For a $C-adsorbate system, got $(length(components)) component counts"))
            all(>=(0), components) || throw(ArgumentError(
                "component counts must be nonnegative, got $components"))
            sum(components) <= length(all_sites) || throw(ArgumentError(
                "component counts request $(sum(components)) occupied sites, " *
                "but the adsorption lattice has only $(length(all_sites))"))
            total_requested = sum(components)
            shuffled = if total_requested == 0
                Int[]
            elseif C == 1 && total_requested == length(all_sites)
                collect(eachindex(all_sites))
            else
                randperm(length(all_sites))[1:total_requested]
            end
            cursor = 1
            for c in 1:C
                n = components[c]
                component_masks[c][shuffled[cursor:(cursor + n - 1)]] .= true
                cursor += n
            end
        else
            length(components) == C || throw(ArgumentError(
                "For a $C-adsorbate system, got $(length(components)) component masks"))
            component_masks = [Vector{Bool}(mask) for mask in components]
        end
        _validate_atomic_components(component_masks, length(all_sites))

        lattice = new{C,G}(lattice_atom, surface, adsorbate_atoms, supercell_dimensions,
                           lattice_constant, lattice_constant_c, periodicity,
                           lattice_positions,
                           num_nearest_neighbors, image_multiplicity,
                           cutoff_radii, neighbors,
                           type_of_sites, all_sites, UInt64(0), component_masks,
                           Int[], 0, false, adsorbate_height,
                           ase_lattice, true, Vector{Bool}[], String[], NaN,
                           UInt64(0))
        lattice.geometry_fingerprint = _atomic_geometry_fingerprint(lattice)
        return sync_ase_lattice!(lattice)
    end
end

"""Content signature for every constructor-time field that defines site geometry."""
function _atomic_geometry_fingerprint(lattice::AtomicLattice)
    state = UInt64(0xcbf29ce484222325)
    mix_byte(h, byte) = (h ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    function mix_word(h, word::UInt64)
        for shift in 0:8:56
            h = mix_byte(h, (word >> shift) & 0xff)
        end
        return h
    end
    function mix_text(h, value)
        bytes = codeunits(string(value))
        h = mix_word(h, UInt64(length(bytes)))
        for byte in bytes
            h = mix_byte(h, byte)
        end
        return h
    end
    mix_int(h, value::Integer) = mix_word(h, reinterpret(UInt64, Int64(value)))
    mix_float(h, value::Float64) = mix_word(h, reinterpret(UInt64, value))

    state = mix_text(state, lattice.lattice_atom)
    state = mix_text(state, lattice.surface)
    for value in lattice.supercell_dimensions
        state = mix_int(state, value)
    end
    state = mix_float(state, lattice.lattice_constant)
    if lattice.lattice_constant_c === nothing
        state = mix_byte(state, 0)
    else
        state = mix_byte(state, 1)
        state = mix_float(state, lattice.lattice_constant_c)
    end
    for value in lattice.periodicity
        state = mix_byte(state, value)
    end
    state = mix_int(state, lattice.num_nearest_neighbors)
    state = mix_byte(state, lattice.image_multiplicity)
    state = mix_int(state, size(lattice.lattice_positions, 1))
    state = mix_int(state, size(lattice.lattice_positions, 2))
    for value in lattice.lattice_positions
        state = mix_float(state, value)
    end
    state = mix_int(state, length(lattice.cutoff_radii))
    for value in lattice.cutoff_radii
        state = mix_float(state, value)
    end
    state = mix_int(state, length(lattice.neighbors))
    for site_shells in lattice.neighbors
        state = mix_int(state, length(site_shells))
        for shell in site_shells
            state = mix_int(state, length(shell))
            for site in shell
                state = mix_int(state, site)
            end
        end
    end
    state = mix_int(state, length(lattice.type_of_sites))
    for site_type in lattice.type_of_sites
        state = mix_text(state, site_type)
    end
    state = mix_int(state, length(lattice.all_sites))
    for (x, y) in lattice.all_sites
        state = mix_float(state, x)
        state = mix_float(state, y)
    end
    return state
end

"""In-process signature of the complete derived ASE frame."""
function _atomic_ase_cache_fingerprint(slab)
    cell = pyconvert(Matrix{Float64}, slab.get_cell())
    positions = pyconvert(Matrix{Float64}, slab.get_positions())
    pbc = pyconvert(Vector{Bool}, slab.get_pbc())
    tags = pyconvert(Vector{Int}, slab.get_tags())
    symbols = pyconvert(Vector{String}, slab.get_chemical_symbols())
    adsorption_cell = try
        pyconvert(Matrix{Float64}, slab.info["adsorbate_info"]["cell"])
    catch
        zeros(Float64, 0, 0)
    end
    return hash((size(cell), Tuple(cell), Tuple(pbc), Tuple(tags),
                 Tuple(symbols), size(positions), Tuple(positions),
                 size(adsorption_cell), Tuple(adsorption_cell)))
end

"""Reject direct mutation of an AtomicLattice's derived Python frame."""
function _validate_atomic_ase_cache(lattice::AtomicLattice)
    if lattice.ase_cache_fingerprint != 0 &&
       _atomic_ase_cache_fingerprint(lattice.ase_lattice) !=
           lattice.ase_cache_fingerprint
        throw(ArgumentError(
            "AtomicLattice.ase_lattice is a derived cache and may not be " *
            "mutated directly; update the lattice state through its Julia API"))
    end
    return nothing
end

function _ensure_atomic_reflection_row!(lattice::AtomicLattice, pivot::Int)
    _validate_atomic_components(lattice)
    _validate_atomic_ase_cache(lattice)
    lattice.reflection_row_built && lattice.reflection_pivot == pivot &&
        return lattice.reflection_row
    cell = pyconvert(Matrix{Float64}, lattice.ase_lattice.get_cell())
    tol = isempty(lattice.cutoff_radii) ?
        maximum(norm(view(cell, k, 1:2)) for k in 1:2) * 1e-8 :
        first(lattice.cutoff_radii) * 1e-6
    lattice.reflection_row = _build_atomic_reflection_row(
        lattice.all_sites, cell, lattice.periodicity, tol, pivot)
    lattice.reflection_pivot = pivot
    lattice.reflection_row_built = true
    return lattice.reflection_row
end

"""
    deepcopy(lattice::AtomicLattice)

Copy an `AtomicLattice`, including an independent Python copy of its cached ASE
frame. PythonCall's `Py` wrapper is otherwise copied without cloning the Python
object it refers to, which lets synchronization or an MLIP evaluation through
one Julia copy mutate another copy's cache.

The component masks are the source of truth. A dirty source produces a
dirty copy whose independent ASE cache is rebuilt on its next
`sync_ase_lattice!` call.
"""
function Base.deepcopy_internal(lattice::AtomicLattice, stackdict::IdDict)
    haskey(stackdict, lattice) && return stackdict[lattice]

    copied = invoke(Base.deepcopy_internal, Tuple{Any, IdDict}, lattice, stackdict)
    copied.ase_lattice = _PY_COPY.deepcopy(lattice.ase_lattice)
    return copied
end

function _validate_atomic_components(
        components::AbstractVector{<:AbstractVector{Bool}}, n_sites::Int)
    all(length(mask) == n_sites for mask in components) || throw(DimensionMismatch(
        "every AtomicLattice component mask must have $n_sites entries"))
    for site in 1:n_sites
        sum(mask[site] for mask in components) <= 1 || throw(ArgumentError(
            "AtomicLattice component masks overlap at site $site"))
    end
    return nothing
end

function _validate_atomic_components(lattice::AtomicLattice{C}) where C
    _atomic_geometry_fingerprint(lattice) == lattice.geometry_fingerprint ||
        throw(ArgumentError(
            "AtomicLattice geometry is fixed at construction; construct a new " *
            "lattice instead of mutating geometry-defining fields"))
    length(lattice.components) == C || throw(DimensionMismatch(
        "AtomicLattice{$C} must contain exactly $C component masks, got " *
        "$(length(lattice.components))"))
    length(lattice.adsorbate_atoms) == C || throw(DimensionMismatch(
        "AtomicLattice{$C} must contain exactly $C adsorbate symbols, got " *
        "$(length(lattice.adsorbate_atoms))"))
    all(!isempty, lattice.adsorbate_atoms) || throw(ArgumentError(
        "AtomicLattice adsorbate symbols must be non-empty"))
    isfinite(lattice.adsorbate_height) || throw(ArgumentError(
        "AtomicLattice adsorbate_height must be finite"))
    return _validate_atomic_components(lattice.components, num_sites(lattice))
end

"""
    coverage(lattice::AtomicLattice)

Fractional coverage derived directly from the component occupation masks.
"""
coverage(lattice::AtomicLattice) =
    sum(sum, lattice.components) / length(lattice.all_sites)

"""
    nn_distance(lattice::AtomicLattice)

Shortest primitive in-plane translation of the selected ASE surface.
The site-finding geometry uses this surface-dependent distance.
"""
function nn_distance(lattice::AtomicLattice)
    _validate_atomic_components(lattice)
    _validate_atomic_ase_cache(lattice)
    if lattice.surface == :fcc211
        cell = pyconvert(Matrix{Float64}, lattice.ase_lattice.get_cell())
        nx, ny, _ = lattice.supercell_dimensions
        # ASE's fcc211 builder groups the first size dimension in triples:
        # a size of `nx` contains `nx / 3` primitive repeats along cell vector 1.
        return min(norm(view(cell, 1, 1:2)) / (nx ÷ 3),
                   norm(view(cell, 2, 1:2)) / ny)
    end
    info = lattice.ase_lattice.info["adsorbate_info"]
    primitive = pyconvert(Matrix{Float64}, info["cell"])
    a = collect(view(primitive, 1, :))
    b = collect(view(primitive, 2, :))
    return minimum((norm(a), norm(b), norm(a - b), norm(a + b)))
end

"""
    sync_ase_lattice!(lattice::AtomicLattice)

Rebuild `ase_lattice`'s adsorbates from `components`, and clear `ase_dirty`.

`ase_lattice` is a cache. Occupancy moves update `components` and set
`ase_dirty`; direct component-mask mutation is detected by comparison with the
last synchronized occupancy snapshot. The ASE frame is only made to agree when
something needs to look at it—an energy evaluation through a Python calculator,
viewing, or writing a trajectory.

Code that needs the current atom list must call this first. Geometry-only cache
readers call `_validate_atomic_ase_cache` so direct Python-side mutation is
rejected without rebuilding adsorbates unnecessarily.

Adsorbates are identified by ASE tag 0, which is what `ase.build.add_adsorbate`
assigns; the substrate keeps the layer tags its ASE surface builder gave it.
Deletion goes in reverse index order because removing an atom renumbers
everything after it.
"""
function sync_ase_lattice!(lattice::AtomicLattice)
    _validate_atomic_components(lattice)
    _validate_atomic_ase_cache(lattice)
    !lattice.ase_dirty && lattice.components == lattice.synced_components &&
        lattice.adsorbate_atoms == lattice.synced_adsorbate_atoms &&
        lattice.adsorbate_height == lattice.synced_adsorbate_height &&
        return lattice

    slab = lattice.ase_lattice
    tags = pyconvert(Vector{Int}, slab.get_tags())
    ads = findall(==(0), tags)
    if !isempty(ads)
        slab.__delitem__(pylist([i - 1 for i in reverse(ads)]))
    end

    for (component, adsorbate) in zip(lattice.components, lattice.adsorbate_atoms)
        for i in findall(component)
            x, y = lattice.all_sites[i]
            ase.build.add_adsorbate(slab, adsorbate;
                                    height=lattice.adsorbate_height, position=(x, y))
        end
    end

    lattice.ase_dirty = false
    lattice.synced_components = copy.(lattice.components)
    lattice.synced_adsorbate_atoms = copy(lattice.adsorbate_atoms)
    lattice.synced_adsorbate_height = lattice.adsorbate_height
    lattice.ase_cache_fingerprint =
        _atomic_ase_cache_fingerprint(lattice.ase_lattice)
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

"""
    _out_of_plane_spacing(lattice_constant, interlayer_spacing)

Resolve the third-axis lattice spacing for the keyword constructors:
`nothing` means isotropic spacing (the in-plane `lattice_constant`); an
explicit value must be finite and strictly positive, the cutoff-ladder
idiom of `compute_neighbors`, and violations throw an `ArgumentError`.
"""
function _out_of_plane_spacing(lattice_constant::Float64, interlayer_spacing::Union{Nothing,Float64})
    interlayer_spacing === nothing && return lattice_constant
    if !(isfinite(interlayer_spacing) && interlayer_spacing > 0.0)
        throw(ArgumentError(
            "interlayer_spacing must be finite and strictly positive, " *
            "got $interlayer_spacing"))
    end
    return interlayer_spacing
end

"""
    _layer_offset(lattice_constant, layer_offset, stacking)

Resolve the in-plane offset between successive layers for the triangular
keyword constructor: `stacking = :aligned` returns `layer_offset` (the zero
offset when `nothing`), and `stacking = :abc` returns the triangle centre
`(a/2, √3·a/6)` of the in-plane lattice, so that each layer sits over the
triangle centres of the one below and the site geometry repeats after
three layers. An explicit `layer_offset` together with `stacking = :abc`,
an unknown `stacking`, or a non-finite offset throws an `ArgumentError`.
No layer-count condition applies: the supercell's third vector is
`supercell_dimensions[3]` times the skewed `a₃`, a lattice vector of the
stacked crystal for every layer count, so every `supercell_dimensions[3]`
is a valid periodic supercell with the same site geometry (see the
[`MLattice`](@ref) docstring).
"""
function _layer_offset(lattice_constant::Float64,
                       layer_offset::Union{Nothing,NTuple{2,Float64}},
                       stacking::Symbol)
    if stacking == :aligned
        offset = layer_offset === nothing ? (0.0, 0.0) : layer_offset
    elseif stacking == :abc
        if layer_offset !== nothing
            throw(ArgumentError(
                "stacking=:abc fixes the layer offset to (a/2, √3·a/6); pass " *
                "either stacking or layer_offset, not both, got " *
                "layer_offset=$layer_offset"))
        end
        offset = (lattice_constant / 2, sqrt(3) * lattice_constant / 6)
    else
        throw(ArgumentError("stacking must be :aligned or :abc, got :$stacking"))
    end
    if !(isfinite(offset[1]) && isfinite(offset[2]))
        throw(ArgumentError("layer_offset must be finite, got $offset"))
    end
    return offset
end

function MLattice{C,SquareLattice}(; lattice_constant::Float64=1.0,
                                    interlayer_spacing::Union{Nothing,Float64}=nothing,
                                    basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0)],
                                    supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 4, 1),
                                    periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                    cutoff_radii::Vector{Float64}=[1.1, 1.5],
                                    components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                    adsorptions::Union{Vector{Int},Symbol}=:full,
                                    image_multiplicity::Bool=false,
                                ) where C

    c = _out_of_plane_spacing(lattice_constant, interlayer_spacing)
    lattice_vectors = [lattice_constant 0.0 0.0; 0.0 lattice_constant 0.0; 0.0 0.0 c]
    lattice_comp, lattice_adsorptions = mlattice_setup(C, basis, supercell_dimensions, components, adsorptions)

    return MLattice{C,SquareLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, cutoff_radii, lattice_comp, lattice_adsorptions;
                                     image_multiplicity=image_multiplicity)
end

function MLattice{C,TriangularLattice}(; lattice_constant::Float64=1.0,
                                        interlayer_spacing::Union{Nothing,Float64}=nothing,
                                        basis::Vector{Tuple{Float64,Float64,Float64}}=[(0.0, 0.0, 0.0),(1/2, sqrt(3)/2, 0.0)],
                                        supercell_dimensions::Tuple{Int64,Int64,Int64}=(4, 2, 1),
                                        periodicity::Tuple{Bool,Bool,Bool}=(true, true, false),
                                        # The triangular second-neighbor distance is √3 ≈ 1.732, so the
                                        # square-lattice default [1.1, 1.5] would leave shell 2 empty.
                                        cutoff_radii::Vector{Float64}=[1.1, 1.8],
                                        components::Union{Vector{Vector{Int64}},Vector{Vector{Bool}},Symbol}=:equal,
                                        adsorptions::Union{Vector{Int},Symbol}=:full,
                                        image_multiplicity::Bool=false,
                                        layer_offset::Union{Nothing,NTuple{2,Float64}}=nothing,
                                        stacking::Symbol=:aligned,
                                    ) where C

    c = _out_of_plane_spacing(lattice_constant, interlayer_spacing)
    dx, dy = _layer_offset(lattice_constant, layer_offset, stacking)
    lattice_vectors = [lattice_constant 0.0 dx; 0.0 sqrt(3)*lattice_constant dy; 0.0 0.0 c]
    lattice_comp, lattice_adsorptions = mlattice_setup(C, basis, supercell_dimensions, components, adsorptions)

    return MLattice{C,TriangularLattice}(lattice_vectors, basis, supercell_dimensions, periodicity, cutoff_radii, lattice_comp, lattice_adsorptions;
                                         image_multiplicity=image_multiplicity)

end





const SLattice{G} = MLattice{1,G} # alias for single-component lattices

const GLattice{C} = MLattice{C,GenericLattice} # alias for generic lattices

num_lattice_components(lattice::MLattice{C,G}) where {C,G} = C

"""
    num_lattice_components(lattice::AtomicLattice{C,G}) where {C,G}

Number of adsorbate species on an `AtomicLattice`, i.e. its first type parameter.
`LatticeWalker` uses this value as its component-count type parameter.
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

occupied_site_count(lattice::AtomicLattice) = sum.(lattice.components)

"""Fractional in-plane coordinates of AtomicLattice adsorption sites."""
function _atomic_site_fractions(lattice::AtomicLattice)
    _validate_atomic_components(lattice)
    _validate_atomic_ase_cache(lattice)
    cell = pyconvert(Matrix{Float64}, lattice.ase_lattice.get_cell())
    inplane = [cell[1, 1] cell[2, 1]; cell[1, 2] cell[2, 2]]
    reciprocal = inv(inplane)
    return [mod.(reciprocal * [x, y], 1.0) for (x, y) in lattice.all_sites]
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

function order_parameter_c2x2(lattice::AtomicLattice{1,SquareLattice})
    d1, d2, _ = lattice.supercell_dimensions
    num_sites(lattice) == d1 * d2 || throw(ArgumentError(
        "order_parameter_c2x2 for AtomicLattice requires one translational " *
        "adsorption-site orbit ($((d1 * d2)) sites), got $(num_sites(lattice))"))
    if isodd(d1) || isodd(d2)
        throw(ArgumentError("order_parameter_c2x2 requires even in-plane " *
            "supercell dimensions, got ($d1, $d2)"))
    end
    return bragg_amplitude(lattice, d1 ÷ 2, d2 ÷ 2)
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
default `supercell_dimensions = (4, 2, 1)` is *not* commensurate. For
multilayer cells use the layer-indexed method
`order_parameter_sqrt3(lattice, layer)`.
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
    _check_triangular_basis(caller::Symbol, lattice::MLattice{1,TriangularLattice})

Shared guard for the layer-indexed triangular observables: throw an
`ArgumentError`, naming the public entry point `caller`, unless the lattice
carries the standard two-site centered-rectangular basis
`[(0, 0, 0), (a/2, √3·a/2, 0)]` consistent with its lattice vectors, the
convention under which the in-plane site decoders of the single-layer
kernels apply.
"""
function _check_triangular_basis(caller::Symbol, lattice::MLattice{1,TriangularLattice})
    if length(lattice.basis) != 2
        throw(ArgumentError("$caller requires the two-site " *
            "centered-rectangular triangular basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    ax, ay = lattice.lattice_vectors[1, 1], lattice.lattice_vectors[2, 2]
    b1, b2 = lattice.basis
    if !(isapprox(b1[1], 0.0, atol=1e-9) && isapprox(b1[2], 0.0, atol=1e-9) &&
         isapprox(b2[1], ax / 2, atol=1e-9) && isapprox(b2[2], ay / 2, atol=1e-9))
        throw(ArgumentError("$caller requires the standard " *
            "triangular basis [(0, 0, 0), (a/2, √3·a/2, 0)] consistent with " *
            "the lattice vectors, got $(lattice.basis)"))
    end
    return nothing
end

"""
    _check_layer_index(caller::Symbol, d3::Int, layer::Int)

Shared guard for the layer-indexed observables: throw an `ArgumentError`,
naming the public entry point `caller`, unless `1 ≤ layer ≤ d3`, in the
style of [`layer_coverage`](@ref).
"""
function _check_layer_index(caller::Symbol, d3::Int, layer::Int)
    if !(1 <= layer <= d3)
        throw(ArgumentError("$caller requires 1 <= layer <= $d3 " *
            "(= supercell_dimensions[3]), got $layer"))
    end
    return nothing
end

"""
    order_parameter_sqrt3(lattice::MLattice{1,TriangularLattice}, layer::Int) -> Float64

Layer-resolved form of [`order_parameter_sqrt3`](@ref) for multilayer
triangular cells: the same three-sublattice modulus evaluated on the
contiguous block of `B = 2·d₁·d₂` sites of layer `layer` (see
[`site_layers`](@ref)) and normalized by `B`, so a perfect √3×√3
arrangement in that layer gives `1/3` whatever the other layers hold.
Because `lattice_positions` orders basis innermost, dimension 1 next and
dimension 3 outermost, a layer block carries the site ordering of the
two-dimensional cell and the tripartition arithmetic of the single-layer
kernel applies verbatim, for aligned and for offset-stacked layers alike
(the `layer_offset`/`stacking` keywords of the constructor): the sublattice
labels are read from the block index, not from the Cartesian positions. At
`d₃ == 1`, `layer = 1` returns exactly the zero-argument value (same loop,
same reduction).

Usable per configuration through the `observables` keyword of the
nested-sampling loops, one callback per layer:

    observables = [Symbol(:psi, k) => (cfg -> order_parameter_sqrt3(cfg, k)) for k in 1:d₃]

Requires the standard two-site basis and `supercell_dimensions[1]`
divisible by 3, as the zero-argument form, plus `1 ≤ layer ≤ d₃`; the
two-dimensionality guard of the zero-argument form does not apply.
Violations throw an `ArgumentError`.
"""
function order_parameter_sqrt3(lattice::MLattice{1,TriangularLattice}, layer::Int)
    _check_triangular_basis(:order_parameter_sqrt3, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    _check_layer_index(:order_parameter_sqrt3, d3, layer)
    if d1 % 3 != 0
        throw(ArgumentError("order_parameter_sqrt3 requires " *
            "supercell_dimensions[1] divisible by 3 so the three √3×√3 " *
            "sublattices close on the periodic cell, got $d1"))
    end
    B = 2 * d1 * d2
    occ = lattice.components[1]
    s0 = (layer - 1) * B
    # The single-layer sublattice decoder applied to the block-local index
    # t = s − s0 (dimension 3 outermost, so the block is one layer)
    n = zeros(Int, 3)
    for t in 1:B
        occ[s0 + t] || continue
        b = (t - 1) % 2
        ci = ((t - 1) ÷ 2) % d1
        c = b == 0 ? ci % 3 : (ci + 2) % 3
        n[c+1] += 1
    end
    s2 = n[1]^2 + n[2]^2 + n[3]^2 - n[1] * n[2] - n[2] * n[3] - n[3] * n[1]
    return sqrt(Float64(s2)) / B
end

"""
    bragg_amplitude(lattice::MLattice{1,SquareLattice}, m::Int, n::Int) -> Float64

Normalized Bragg-peak amplitude of the occupation pattern on a
single-component square lattice:

    |ρ(k)| = |Σ_{occupied sites} e^{ik·r}| / M,   k = 2π·(m/d₁, n/d₂),

in units of the inverse lattice constant, where `(d₁, d₂)` are the in-plane
supercell dimensions and `M = d₁·d₂` is the number of sites. This is the
LEED-intensity convention of Zhang, Blum & Reuter
[PRB **75**, 235406 (2007)]; [`order_parameter_c2x2`](@ref) is the
k = (π, π) member of the family, and
`bragg_amplitude(lat, d1 ÷ 2, d2 ÷ 2) == order_parameter_c2x2(lat)` on
even-dimension cells. Integer indices on the supercell reciprocal grid make
every representable k commensurate by construction; indices wrap modulo the
dimensions, so negative and out-of-range values are valid:
`bragg_amplitude(lat, m, n) == bragg_amplitude(lat, m + d1, n)`. The empty
and the full lattice give `0` for any `(m, n)` not both `≡ 0` modulo the
dimensions; a single particle gives `1/M`.

Like the named order parameters, |ρ(k)| is not a function of `(E, N)` and
must be evaluated on configurations (e.g. per culled walker, via the
`observables` keyword of the nested-sampling loops); see
[`order_parameter_stripe`](@ref) for the composed-observable recipes.

Requires a single-site basis and a strictly two-dimensional supercell
(`supercell_dimensions[3] == 1`); violations throw an `ArgumentError`.
"""
function bragg_amplitude(lattice::MLattice{1,SquareLattice}, m::Int, n::Int)
    if length(lattice.basis) != 1
        throw(ArgumentError("bragg_amplitude requires a single-site " *
            "basis, got $(length(lattice.basis)) basis sites"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("bragg_amplitude requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    # Phase tables over mod-reduced cispi arguments: half-turn phases are
    # exactly ±1.0 and quarter-turn phases exactly ±1.0/±1.0im, so the c2x2
    # identity and the period-4 values hold to full precision
    col = [cispi(2 * mod(m * x, d1) / d1) for x in 0:d1-1]
    row = [cispi(2 * mod(n * y, d2) / d2) for y in 0:d2-1]
    occ = lattice.components[1]
    z = 0.0 + 0.0im
    for s in eachindex(occ)
        if occ[s]
            # lattice_positions ordering: basis innermost, dimension 1 fastest
            i0 = (s - 1) % d1
            j0 = ((s - 1) ÷ d1) % d2
            z += col[i0+1] * row[j0+1]
        end
    end
    return abs(z) / (d1 * d2)
end

function bragg_amplitude(lattice::AtomicLattice{1}, m::Int, n::Int)
    lattice.periodicity[1] && lattice.periodicity[2] || throw(ArgumentError(
        "bragg_amplitude for AtomicLattice requires periodic x and y boundaries"))
    fractions = _atomic_site_fractions(lattice)
    z = 0.0 + 0.0im
    for site in occupied_indices(lattice)
        f = fractions[site]
        z += cispi(2 * (m * f[1] + n * f[2]))
    end
    return abs(z) / num_sites(lattice)
end

"""
    order_parameter_sqrt3(lattice::AtomicLattice{1,TriangularLattice}) -> Float64

Three-sublattice order parameter for a one-site triangular adsorption lattice,
evaluated at a K point of the primitive ASE surface cell.  A perfect
(sqrt(3) x sqrt(3))R30-degree overlayer at coverage 1/3 returns `1/3`.

The selected adsorption-site family must form exactly one translational orbit
(`num_sites(lattice) == d1*d2`), and both in-plane dimensions must be divisible
by three so the ordered state closes across the periodic boundaries.
"""
function order_parameter_sqrt3(lattice::AtomicLattice{1,TriangularLattice})
    d1, d2, _ = lattice.supercell_dimensions
    num_sites(lattice) == d1 * d2 || throw(ArgumentError(
        "order_parameter_sqrt3 for AtomicLattice requires one translational " *
        "adsorption-site orbit ($((d1 * d2)) sites), got $(num_sites(lattice))"))
    if d1 % 3 != 0 || d2 % 3 != 0
        throw(ArgumentError("order_parameter_sqrt3 for AtomicLattice requires " *
            "both in-plane dimensions divisible by 3, got ($d1, $d2)"))
    end
    return bragg_amplitude(lattice, d1 ÷ 3, -(d2 ÷ 3))
end

"""
    order_parameter_p2x2(lattice::AtomicLattice{1,TriangularLattice}) -> Float64

Orientation-independent p(2x2) order parameter for a one-site triangular
adsorption lattice.  It is the quadrature sum of the three primitive-cell M
points; a perfect p(2x2) overlayer at coverage 1/4 returns `sqrt(3)/4`.
"""
function order_parameter_p2x2(lattice::AtomicLattice{1,TriangularLattice})
    d1, d2, _ = lattice.supercell_dimensions
    num_sites(lattice) == d1 * d2 || throw(ArgumentError(
        "order_parameter_p2x2 for AtomicLattice requires one translational " *
        "adsorption-site orbit ($((d1 * d2)) sites), got $(num_sites(lattice))"))
    if isodd(d1) || isodd(d2)
        throw(ArgumentError("order_parameter_p2x2 for AtomicLattice requires " *
            "even in-plane dimensions, got ($d1, $d2)"))
    end
    return sqrt(bragg_amplitude(lattice, d1 ÷ 2, 0)^2 +
                bragg_amplitude(lattice, 0, d2 ÷ 2)^2 +
                bragg_amplitude(lattice, d1 ÷ 2, d2 ÷ 2)^2)
end

"""
    order_parameter_stripe(lattice::MLattice{1,SquareLattice}; period::Int = 2) -> Float64

Orientation-degenerate axial stripe order parameter on a single-component
square lattice: the quadrature sum of the two axial Bragg amplitudes at
wavevector 2π/P, with P = `period` in lattice constants,

    Ψ = √(|ρ(2π/P, 0)|² + |ρ(0, 2π/P)|²),

i.e. `sqrt(bragg_amplitude(lat, d1 ÷ period, 0)^2 +
bragg_amplitude(lat, 0, d2 ÷ period)^2)` with the normalized amplitudes of
[`bragg_amplitude`](@ref). A perfect single-orientation period-P stripe
with the half-period filled attains `1/(P·sin(π/P))`: each of the `d₂` rows
contributes `d₁/P` copies of the geometric phasor sum
`Σ_{x=0}^{P/2−1} e^{2πix/P}`, of modulus `1/sin(π/P)`, and dividing by
`M = d₁·d₂` leaves `1/(P·sin(π/P))`. That gives `1/2` at P = 2 (matching
the c(2×2) = 1/2 convention at half filling) and `√2/4` at P = 4; the
perpendicular component vanishes on a perfect stripe, so the quadrature
maximum equals the single-orientation value. `period = 2` detects the (2×1)
row and column stripes — the superantiferromagnetic phase of the square
lattice with nearest- and next-nearest-neighbor couplings [Binder & Landau,
PRB **21**, 1941 (1980)] — `period = 4` the axial period-4 phases that
appear with third-neighbor couplings [Landau & Binder, PRB **31**, 5946
(1985)], and `period = 2h` the width-h stripes of dipolar-frustrated
ferromagnets [MacIsaac, Whitehead, Robinson & De'Bell, PRB **51**, 16033
(1995)]. The perfect checkerboard and the empty and the full lattice
give `0`.

Ψ is not a function of `(E, N)` — degenerate energy levels contain ordered
and disordered configurations alike — so it must be evaluated on
configurations (e.g. per culled walker, via the `observables` keyword of
the nested-sampling loops). Higher moments for Binder-cumulant analysis are
composed caller-side, with no library change:

    observables = [:stripe  => order_parameter_stripe,
                   :stripe2 => cfg -> order_parameter_stripe(cfg)^2,
                   :stripe4 => cfg -> order_parameter_stripe(cfg)^4]

Diagonal period-4 order is composed the same way, from the k = (π/2, ±π/2)
quadrature

    cfg -> sqrt(bragg_amplitude(cfg, d1 ÷ 4, d2 ÷ 4)^2 +
                bragg_amplitude(cfg, d1 ÷ 4, -(d2 ÷ 4))^2)

which equals `√2/4` on every member of the k = (π/2, ±π/2) ground manifold
of the nearest-neighbor-attraction plus isotropic third-neighbor-repulsion
model — p(2×2) blocks and diagonal period-4 stripes alike — while
`order_parameter_stripe(period=4)` is exactly `0` there; the two
observables separate the axial and diagonal period-4 phases.

Requires a single-site basis, a strictly two-dimensional supercell
(`supercell_dimensions[3] == 1`), `period >= 2`, and `period` dividing both
in-plane dimensions (both orientations are evaluated, and an incommensurate
orientation leaks a spurious amplitude); violations throw an
`ArgumentError`.
"""
function order_parameter_stripe(lattice::MLattice{1,SquareLattice}; period::Int=2)
    if length(lattice.basis) != 1
        throw(ArgumentError("order_parameter_stripe requires a single-site " *
            "basis, got $(length(lattice.basis)) basis sites"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("order_parameter_stripe requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    if period < 2
        throw(ArgumentError("order_parameter_stripe requires period >= 2, " *
            "got $period"))
    end
    if d1 % period != 0 || d2 % period != 0
        bad = d1 % period != 0 ? d1 : d2
        throw(ArgumentError("order_parameter_stripe requires the period to " *
            "divide both in-plane supercell dimensions, since both stripe " *
            "orientations are evaluated and an incommensurate orientation " *
            "leaks a spurious amplitude; period $period does not divide $bad"))
    end
    return sqrt(bragg_amplitude(lattice, d1 ÷ period, 0)^2 +
                bragg_amplitude(lattice, 0, d2 ÷ period)^2)
end

function order_parameter_stripe(lattice::AtomicLattice{1,SquareLattice}; period::Int=2)
    d1, d2, _ = lattice.supercell_dimensions
    num_sites(lattice) == d1 * d2 || throw(ArgumentError(
        "order_parameter_stripe for AtomicLattice requires one translational " *
        "adsorption-site orbit ($((d1 * d2)) sites), got $(num_sites(lattice))"))
    period >= 2 || throw(ArgumentError(
        "order_parameter_stripe requires period >= 2, got $period"))
    if d1 % period != 0 || d2 % period != 0
        bad = d1 % period != 0 ? d1 : d2
        throw(ArgumentError("order_parameter_stripe requires period $period " *
            "to divide both in-plane dimensions; it does not divide $bad"))
    end
    return sqrt(bragg_amplitude(lattice, d1 ÷ period, 0)^2 +
                bragg_amplitude(lattice, 0, d2 ÷ period)^2)
end

"""
    bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int) -> Float64

Normalized Bragg-peak amplitude of the occupation pattern on a
single-component triangular lattice:

    |ρ(k)| = |Σ_{occupied sites} e^{ik·r}| / M,   k = 2π·(m/d₁, n/(√3·d₂)),

in units of the inverse lattice constant, where `(d₁, d₂)` are the in-plane
supercell dimensions of the conventional (centered-rectangular) cell and
`M = 2·d₁·d₂` is the number of sites; the sum runs over both basis sites,
so the two-site basis phases are included. Integer indices on the
conventional-cell reciprocal grid make every representable k commensurate
by construction, in the same LEED-intensity convention as the square-lattice
method [Zhang, Blum & Reuter, PRB **75**, 235406 (2007)].

Because the basis sites live on the half-integer grid of the conventional
cell, amplitudes are periodic under the *primitive* reciprocal lattice
rather than index-wise:
`bragg_amplitude(lat, m, n) == bragg_amplitude(lat, m + d1, n - d2) ==
bragg_amplitude(lat, m, n + 2*d2)` (the two index shifts add the primitive
reciprocal vectors b₁ = 2π·(1, −1/√3) and b₂ = 2π·(0, 2/√3)), and both
identities hold exactly. The three M points of the triangular Brillouin
zone sit at `(m, n) = (0, d2)` (representable on any cell) and
`(d1 ÷ 2, ∓(d2 ÷ 2))` (integer indices for even in-plane dimensions);
M-point phases are exactly ±1.0 (u and v share the parity of the basis
index, so the quarter-turn table entries always combine to half turns), so
the documented values of [`order_parameter_p2x2`](@ref) hold to full
precision. The empty lattice gives exactly `0` at every k; the full
lattice gives `0` for any k that is not a reciprocal-lattice vector of the
triangular lattice, exactly at the three M points and to floating-point
roundoff at generic indices (where the phasor components are irrational).
A single particle gives `1/M`.

Like the named order parameters, |ρ(k)| is not a function of `(E, N)` and
must be evaluated on configurations (e.g. per culled walker, via the
`observables` keyword of the nested-sampling loops); see
[`order_parameter_p2x2`](@ref) for the composed-observable recipes.

Requires the standard two-site centered-rectangular triangular basis
`[(0, 0, 0), (a/2, √3·a/2, 0)]` consistent with the lattice vectors and a
strictly two-dimensional supercell (`supercell_dimensions[3] == 1`);
violations throw an `ArgumentError`. For multilayer cells use the
layer-indexed method `bragg_amplitude(lattice, m, n, layer)` or
[`bragg_amplitude_layers`](@ref).
"""
function bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int)
    if length(lattice.basis) != 2
        throw(ArgumentError("bragg_amplitude requires the two-site " *
            "centered-rectangular triangular basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    ax, ay = lattice.lattice_vectors[1, 1], lattice.lattice_vectors[2, 2]
    b1, b2 = lattice.basis
    if !(isapprox(b1[1], 0.0, atol=1e-9) && isapprox(b1[2], 0.0, atol=1e-9) &&
         isapprox(b2[1], ax / 2, atol=1e-9) && isapprox(b2[2], ay / 2, atol=1e-9))
        throw(ArgumentError("bragg_amplitude requires the standard " *
            "triangular basis [(0, 0, 0), (a/2, √3·a/2, 0)] consistent with " *
            "the lattice vectors, got $(lattice.basis)"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("bragg_amplitude requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    # Site (b, ci, cj) sits at r = ((ci − 1) + b/2, √3·((cj − 1) + b/2))·a,
    # so on the half-integer grid u = 2(ci − 1) + b, v = 2(cj − 1) + b the
    # phase is k·r = π·(m·u·d₂ + n·v·d₁)/(d₁·d₂). The numerator is reduced
    # as one integer modulo 2·d₁·d₂ before the single cispi call per
    # occupied site: half- and quarter-turn phases come out exactly
    # ±1.0/±1.0im, every combined M-point phase is exactly ±1.0, and the
    # joint reduction makes the b₁/b₂ wrapping identities hold bit-exactly,
    # since both shifts move the numerator by even multiples of d₁·d₂
    # (u and v share the parity of the basis index).
    twoM = 2 * d1 * d2
    mr = mod(m, 2 * d1)
    nr = mod(n, 2 * d2)
    colnum = [mod(mr * u * d2, twoM) for u in 0:2*d1-1]
    rownum = [mod(nr * v * d1, twoM) for v in 0:2*d2-1]
    occ = lattice.components[1]
    z = 0.0 + 0.0im
    for s in eachindex(occ)
        if occ[s]
            # lattice_positions ordering: basis innermost, dimension 1 fastest
            b = (s - 1) % 2
            ci0 = ((s - 1) ÷ 2) % d1
            cj0 = (s - 1) ÷ (2 * d1)
            num = colnum[2*ci0+b+1] + rownum[2*cj0+b+1]
            num >= twoM && (num -= twoM)
            z += cispi(num / (d1 * d2))
        end
    end
    return abs(z) / twoM
end

"""
    bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int, layer::Int) -> Float64

Layer-resolved form of the triangular [`bragg_amplitude`](@ref) for
multilayer cells: `|ρ(k)|` of the occupation pattern in layer `layer`'s
contiguous block of `B = 2·d₁·d₂` sites (see [`site_layers`](@ref)),
normalized by `B`, at the same in-plane `k = 2π·(m/d₁, n/(√3·d₂))` of the
conventional cell. The phase is read from the block-local site index, so
the value is that of the layer's pattern in its own frame; a rigid in-plane
offset of the whole layer (the `layer_offset`/`stacking` keywords of the
constructor) multiplies ρ by a global phase that the modulus drops, so
aligned and offset-stacked cells give the same value. The relative phases
between layers, which carry the stacking sequence of an in-plane order, are
available from [`bragg_amplitude_layers`](@ref). At `d₃ == 1`, `layer = 1`
returns exactly the zero-argument value (same loop, same reduction).

Requires the standard two-site basis and `1 ≤ layer ≤ d₃`; violations throw
an `ArgumentError`.
"""
function bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int, layer::Int)
    _check_triangular_basis(:bragg_amplitude, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    _check_layer_index(:bragg_amplitude, d3, layer)
    colnum, rownum, twoM = _bragg_phase_tables(d1, d2, m, n)
    z = _bragg_block_sum(lattice.components[1], d1, d2, colnum, rownum, twoM,
                         (layer - 1) * twoM)
    return abs(z) / twoM
end

# Own-frame phase tables of the single-layer triangular kernel for the
# in-plane index (m, n): the reduced column and row numerators and the
# modulus 2·d₁·d₂ (see the single-layer `bragg_amplitude` for the derivation)
function _bragg_phase_tables(d1::Int, d2::Int, m::Int, n::Int)
    twoM = 2 * d1 * d2
    mr = mod(m, 2 * d1)
    nr = mod(n, 2 * d2)
    colnum = [mod(mr * u * d2, twoM) for u in 0:2*d1-1]
    rownum = [mod(nr * v * d1, twoM) for v in 0:2*d2-1]
    return colnum, rownum, twoM
end

# Unnormalized own-frame Bragg sum over the block of twoM sites following
# site s0: the single-layer loop and reduction applied to the block-local
# index t = s − s0 (dimension 3 outermost, so the block is one layer)
function _bragg_block_sum(occ::Vector{Bool}, d1::Int, d2::Int,
                          colnum::Vector{Int}, rownum::Vector{Int}, twoM::Int, s0::Int)
    z = 0.0 + 0.0im
    for t in 1:twoM
        if occ[s0 + t]
            b = (t - 1) % 2
            ci0 = ((t - 1) ÷ 2) % d1
            cj0 = (t - 1) ÷ (2 * d1)
            num = colnum[2*ci0+b+1] + rownum[2*cj0+b+1]
            num >= twoM && (num -= twoM)
            z += cispi(num / (d1 * d2))
        end
    end
    return z
end

"""
    order_parameter_p2x2(lattice::MLattice{1,TriangularLattice}) -> Float64

Orientation-degenerate M-point order parameter for p(2×2) and p(2×1)/row
ordering on a single-component triangular lattice: the quadrature sum of
the three normalized M-point Bragg amplitudes,

    Ψ = √(|ρ(M₁)|² + |ρ(M₂)|² + |ρ(M₃)|²),

i.e. `sqrt(bragg_amplitude(lat, d1 ÷ 2, -(d2 ÷ 2))^2 +
bragg_amplitude(lat, 0, d2)^2 + bragg_amplitude(lat, d1 ÷ 2, d2 ÷ 2)^2)`
with the amplitudes of [`bragg_amplitude`](@ref). A perfect p(2×2)
arrangement at coverage 1/4 contributes exactly `1/4` at every M point,
giving `√3/4 ≈ 0.4330`; a perfect single-orientation p(2×1) row phase at
coverage 1/2 puts exactly `1/2` on one M point and `0` on the other two,
giving `1/2`; the perfect (√3×√3)R30° state (K-point order) and the empty
and the full lattice give exactly `0`; disordered configurations give
O(1/√M). Translation and orientation degeneracies are divided out, so all
degenerate ordered states of each phase give the same value.

The M-point pattern separates the two phases where the scalar cannot: the
max-to-quadrature ratio of the three amplitudes is `1/√3` for p(2×2)
(three equal M points) and `1` for a single p(2×1) orientation; compose it
caller-side from `bragg_amplitude` when needed. [`order_parameter_sqrt3`](@ref)
is exactly `0` on both M-point phases and this function is exactly `0` on
the √3×√3 phase, so the two observables are orthogonal discriminators.

Ψ is not a function of `(E, N)` and must be evaluated on configurations
(e.g. per culled walker, via the `observables` keyword of the
nested-sampling loops). Higher moments for Binder-cumulant analysis are
composed caller-side, with no library change:

    observables = [:psi  => order_parameter_p2x2,
                   :psi2 => cfg -> order_parameter_p2x2(cfg)^2,
                   :psi4 => cfg -> order_parameter_p2x2(cfg)^4]

Requires the structural guards of [`bragg_amplitude`](@ref) plus even
in-plane supercell dimensions: the M points sit at half-integer reciprocal
indices, and even circumferences are exactly the wrap-invariance condition
of the p(2×2) sublattice. A cell hosting this order parameter and
[`order_parameter_sqrt3`](@ref) simultaneously needs
`supercell_dimensions[1]` divisible by 6 with an even second dimension;
the shipped default `(4, 2, 1)` satisfies this function's guards but not
the √3×√3 one's. For multilayer cells use the layer-indexed method
`order_parameter_p2x2(lattice, layer)`.
"""
function order_parameter_p2x2(lattice::MLattice{1,TriangularLattice})
    if length(lattice.basis) != 2
        throw(ArgumentError("order_parameter_p2x2 requires the two-site " *
            "centered-rectangular triangular basis, got " *
            "$(length(lattice.basis)) basis sites"))
    end
    ax, ay = lattice.lattice_vectors[1, 1], lattice.lattice_vectors[2, 2]
    b1, b2 = lattice.basis
    if !(isapprox(b1[1], 0.0, atol=1e-9) && isapprox(b1[2], 0.0, atol=1e-9) &&
         isapprox(b2[1], ax / 2, atol=1e-9) && isapprox(b2[2], ay / 2, atol=1e-9))
        throw(ArgumentError("order_parameter_p2x2 requires the standard " *
            "triangular basis [(0, 0, 0), (a/2, √3·a/2, 0)] consistent with " *
            "the lattice vectors, got $(lattice.basis)"))
    end
    d1, d2, d3 = lattice.supercell_dimensions
    if d3 != 1
        throw(ArgumentError("order_parameter_p2x2 requires a two-dimensional " *
            "supercell (supercell_dimensions[3] == 1), got $d3 layers"))
    end
    if d1 % 2 != 0 || d2 % 2 != 0
        bad = d1 % 2 != 0 ? d1 : d2
        throw(ArgumentError("order_parameter_p2x2 requires even in-plane " *
            "supercell dimensions, since the M points sit at half-integer " *
            "reciprocal indices and the p(2×2) sublattice closes on the " *
            "periodic cell only for even circumferences; got $bad"))
    end
    return sqrt(bragg_amplitude(lattice, d1 ÷ 2, -(d2 ÷ 2))^2 +
                bragg_amplitude(lattice, 0, d2)^2 +
                bragg_amplitude(lattice, d1 ÷ 2, d2 ÷ 2)^2)
end

"""
    order_parameter_p2x2(lattice::MLattice{1,TriangularLattice}, layer::Int) -> Float64

Layer-resolved form of [`order_parameter_p2x2`](@ref) for multilayer
triangular cells: the quadrature of the three layer-indexed M-point
amplitudes of [`bragg_amplitude`](@ref) for layer `layer`, so a perfect
p(2×1) row phase in that layer gives `1/2` and a perfect p(2×2) arrangement
`√3/4`, whatever the other layers hold and for aligned and offset-stacked
layers alike. At `d₃ == 1`, `layer = 1` returns exactly the zero-argument
value.

Usable per configuration through the `observables` keyword of the
nested-sampling loops, one callback per layer:

    observables = [Symbol(:psi, k) => (cfg -> order_parameter_p2x2(cfg, k)) for k in 1:d₃]

Requires the standard two-site basis and even in-plane supercell
dimensions, as the zero-argument form, plus `1 ≤ layer ≤ d₃`; violations
throw an `ArgumentError`.
"""
function order_parameter_p2x2(lattice::MLattice{1,TriangularLattice}, layer::Int)
    _check_triangular_basis(:order_parameter_p2x2, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    _check_layer_index(:order_parameter_p2x2, d3, layer)
    if d1 % 2 != 0 || d2 % 2 != 0
        bad = d1 % 2 != 0 ? d1 : d2
        throw(ArgumentError("order_parameter_p2x2 requires even in-plane " *
            "supercell dimensions, since the M points sit at half-integer " *
            "reciprocal indices and the p(2×2) sublattice closes on the " *
            "periodic cell only for even circumferences; got $bad"))
    end
    return sqrt(bragg_amplitude(lattice, d1 ÷ 2, -(d2 ÷ 2), layer)^2 +
                bragg_amplitude(lattice, 0, d2, layer)^2 +
                bragg_amplitude(lattice, d1 ÷ 2, d2 ÷ 2, layer)^2)
end

"""
    bragg_amplitude_layers(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int) -> Vector{ComplexF64}

Complex per-layer Bragg amplitudes of the occupation pattern on a
multilayer single-component triangular lattice, each layer in its own
frame,

    ρ_k = (1/B) Σ_{occupied sites s in layer k} e^{i k·r_s},   k = 1:d₃,

with `B = 2·d₁·d₂` sites per layer, the in-plane wavevector
`k = 2π·(m/d₁, n/(√3·d₂))` in inverse lattice constants (the convention of
the triangular [`bragg_amplitude`](@ref)), and `r_s` the in-plane position
of site `s` relative to its own layer's origin: the phase table of the
layer-indexed [`bragg_amplitude`](@ref), read from the block-local site
index, so `abs.(ρ)` equals that function's values and the relative phases
between layers compare each layer's pattern in the layer's own frame. The
rigid in-plane offset of an offset-stacked cell (the
`layer_offset`/`stacking` keywords of the constructor) does not enter: it
is absorbed in the skewed supercell whose reciprocal lattice
[`stacking_bragg_amplitude`](@ref) samples, whereas a per-layer offset
phase would make the layer transform depend on which layer is called
layer 1 whenever `d₃·k·Δ` is not a whole number of turns. The vector is
translation-covariant: an in-plane lattice translation multiplies every
`ρ_k` by one common phase, and a translation by the third lattice vector
relabels the layers cyclically.

The return value is a `Vector`, not a scalar, so this function is not
directly usable as an `observables` callback (callbacks must return a
`Real`); record a component's modulus caller-side, e.g.
`cfg -> abs(bragg_amplitude_layers(cfg, m, n)[k])`, or use
[`stacking_bragg_amplitude`](@ref).

Requires the standard two-site centered-rectangular triangular basis
`[(0, 0, 0), (a/2, √3·a/2, 0)]` consistent with the lattice vectors;
violations throw an `ArgumentError`. Any number of layers is accepted,
`d₃ == 1` included.
"""
function bragg_amplitude_layers(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int)
    _check_triangular_basis(:bragg_amplitude_layers, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    colnum, rownum, twoM = _bragg_phase_tables(d1, d2, m, n)
    occ = lattice.components[1]
    rho = zeros(ComplexF64, d3)
    for k in 1:d3
        rho[k] = _bragg_block_sum(occ, d1, d2, colnum, rownum, twoM, (k - 1) * twoM) / twoM
    end
    return rho
end

"""
    stacking_bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int, l::Int) -> Float64

Stacking-resolved Bragg amplitude of a multilayer single-component
triangular lattice: the modulus of the discrete Fourier transform, along
the layer index, of the own-frame per-layer amplitudes of
[`bragg_amplitude_layers`](@ref),

    |A_l| = |(1/d₃) Σ_{k=1}^{d₃} ρ_k e^{2πi·l·(k − 1)/d₃}|,   l ∈ 0:d₃−1.

This is the three-dimensional structure factor of the periodic supercell,
`|Σ_{occupied sites} e^{i K·r}| / M` with `M = 2·d₁·d₂·d₃`, at its own
reciprocal-lattice point `K = 2π·A⁻ᵀ·(m, n, l)` (`A` the supercell
matrix, so the skewed third vector of an offset-stacked cell is accounted
for exactly). It is therefore invariant under every lattice translation
of the configuration, in-plane or along the stacking direction, at every
`(m, n, l)` on every cell, and the `l` labels do not depend on the sign or
size of the layer offset. `l = 0` is the amplitude of the sequence in
which every layer carries the same own-frame pattern, and `l ≥ 1` picks
out the sequences whose own-frame phase advances by `−2π·l/d₃` per layer:
on a three-layer cell at a K point, a √3×√3 sublattice held fixed across
layers gives `l = 0` and one that rotates by one sublattice per layer
gives `l = 1` or `l = 2` by its sense; at an M point, rows held fixed give
`l = 0`, while rows shifted by one row per layer (an own-frame phase
advance of π, a two-layer period) give a pure `l = d₃/2` component on an
even layer count and spread over every `l`, `l = 0` included, on an odd
one ((1/6, 1/3, 1/3) on three layers). The components
obey the Parseval identity `Σ_l |A_l|² = (1/d₃) Σ_k |ρ_k|²`, and each lies
in `[0, 1]`. The value is a scalar, so
`cfg -> stacking_bragg_amplitude(cfg, m, n, l)` is directly usable as an
`observables` callback in the nested-sampling loops.

Requires the structural guard of [`bragg_amplitude_layers`](@ref) and
`0 ≤ l ≤ d₃ − 1`; violations throw an `ArgumentError`.
"""
function stacking_bragg_amplitude(lattice::MLattice{1,TriangularLattice}, m::Int, n::Int, l::Int)
    d3 = lattice.supercell_dimensions[3]
    if !(0 <= l <= d3 - 1)
        throw(ArgumentError("stacking_bragg_amplitude requires 0 <= l <= $(d3 - 1) " *
            "(= supercell_dimensions[3] - 1), got $l"))
    end
    rho = bragg_amplitude_layers(lattice, m, n)
    z = 0.0 + 0.0im
    for k in 1:d3
        # Layer phase reduced modulo one turn, so half turns are exactly ±1
        z += rho[k] * cispi(2 * mod(l * (k - 1), d3) / d3)
    end
    return abs(z) / d3
end

"""
    _check_planar_basis(caller::Symbol, lattice::MLattice)

Shared guard for the layer helpers: throw an `ArgumentError`, naming the
public entry point `caller`, unless all basis z-components are equal.
Layers along dimension 3 are geometrically well defined only for a planar
basis.
"""
function _check_planar_basis(caller::Symbol, lattice::MLattice)
    zs = [b[3] for b in lattice.basis]
    if any(z -> z != zs[1], zs)
        throw(ArgumentError("$caller requires a planar basis (all basis " *
            "z-components equal) so layers along dimension 3 are " *
            "geometrically well defined, got basis z-components $zs"))
    end
    return nothing
end

"""
    site_layers(lattice::MLattice) -> Vector{Int}

Layer index of every site along the third supercell dimension: site `s`
belongs to layer `(s − 1) ÷ B + 1` with `B = length(basis) · d₁ · d₂`,
exact because `lattice_positions` orders sites with dimension 3 outermost
(basis innermost, dimension 1 fastest), making each layer one contiguous
block of `B` sites. The result has `num_sites(lattice)` entries with
values in `1:d₃`.

Requires a planar basis (all basis z-components equal), so that "layer" is
geometrically unambiguous; violations throw an `ArgumentError`.
"""
function site_layers(lattice::MLattice)
    _check_planar_basis(:site_layers, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    B = length(lattice.basis) * d1 * d2
    # lattice_positions ordering: dimension 3 outermost, so layers are
    # contiguous blocks of B sites
    return [(s - 1) ÷ B + 1 for s in 1:(B * d3)]
end

"""
    layer_coverage(lattice::MLattice{1,G}, layer::Int) -> Float64

Occupied fraction of one layer of a single-component lattice: the number
of occupied sites in layer `layer`'s contiguous block of
`B = length(basis) · d₁ · d₂` sites (see [`site_layers`](@ref)), divided
by `B`. This is the layer-resolved complement of the in-plane order
parameters for three-dimensional supercells; the per-layer coverages
θ_k are the order parameters of lattice-gas layering transitions. It
returns a scalar `Real`, so `cfg -> layer_coverage(cfg, k)`
is directly usable as an `observables` callback in the nested-sampling
loops: like the shipped order parameters, θ is not a function of `(E, N)`
and must be evaluated per configuration rather than reconstructed from an
energy ledger. At `d₃ == 1` it degenerates to the total coverage `N/M`,
which is legal and intentional, so no two-dimensionality guard applies.

Requires a planar basis (all basis z-components equal) and
`1 ≤ layer ≤ d₃`; violations throw an `ArgumentError`.
"""
function layer_coverage(lattice::MLattice{1,G}, layer::Int) where G
    _check_planar_basis(:layer_coverage, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    if !(1 <= layer <= d3)
        throw(ArgumentError("layer_coverage requires 1 <= layer <= $d3 " *
            "(= supercell_dimensions[3]), got $layer"))
    end
    B = length(lattice.basis) * d1 * d2
    occ = lattice.components[1]
    n = 0
    for s in ((layer - 1) * B + 1):(layer * B)
        n += occ[s] ? 1 : 0
    end
    return n / B
end

"""
    occupancy_profile(lattice::MLattice{1,G}) -> Vector{Float64}

Vector of layer coverages, `[layer_coverage(lattice, k) for k in 1:d₃]`,
computed in one pass. Every layer holds `B = length(basis) · d₁ · d₂`
sites, so `mean(occupancy_profile(lat)) == N/M` exactly (`N` occupied
sites, `M` total sites).

The return value is a `Vector`, not a scalar, so `occupancy_profile` is
*not* directly usable as an `observables` callback (callbacks must return
a `Real`); record per-layer coverages by composing scalar callbacks
caller-side, with no library change:

    observables = [Symbol(:theta, k) => (cfg -> layer_coverage(cfg, k)) for k in 1:d₃]

after which ⟨θ₁⟩ … ⟨θ_d₃⟩ come from `observable_cols=[:theta1, :theta2, …]`
in the grand-canonical stats functions.

Requires a planar basis (all basis z-components equal); violations throw
an `ArgumentError`.
"""
function occupancy_profile(lattice::MLattice{1,G}) where G
    _check_planar_basis(:occupancy_profile, lattice)
    d1, d2, d3 = lattice.supercell_dimensions
    B = length(lattice.basis) * d1 * d2
    occ = lattice.components[1]
    counts = zeros(Int, d3)
    for s in eachindex(occ)
        if occ[s]
            counts[(s - 1) ÷ B + 1] += 1
        end
    end
    return counts ./ B
end

"""
    layer_field(lattice::MLattice, per_layer::AbstractVector) -> Vector

Broadcast a per-layer value over every site of its layer: the result has
`num_sites(lattice)` entries, with `per_layer[k]` at every site of layer
`k` (see [`site_layers`](@ref)). The element type of `per_layer` is
preserved, so a `Unitful` profile feeds a `SiteFieldLatticeHamiltonian`
directly:

    field = layer_field(lat, [-0.27, -0.03375, -0.01] .* u"eV")
    h = SiteFieldLatticeHamiltonian(GenericLatticeHamiltonian(0.0, [-0.01], u"eV"), field)

A height profile is physically meaningful only when dimension 3 is
non-periodic (`periodicity[3] == false`); this is not checked.

Requires a planar basis (all basis z-components equal) and
`length(per_layer) == d₃` (= `supercell_dimensions[3]`); violations throw
an `ArgumentError`.
"""
function layer_field(lattice::MLattice, per_layer::AbstractVector)
    _check_planar_basis(:layer_field, lattice)
    d3 = lattice.supercell_dimensions[3]
    if length(per_layer) != d3
        throw(ArgumentError("layer_field requires one value per layer " *
            "(length(per_layer) == supercell_dimensions[3] == $d3), got " *
            "$(length(per_layer)) values"))
    end
    return per_layer[site_layers(lattice)]
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
function enumerate_motif_embeddings(lattice::AbstractLattice, distances::AbstractVector{<:Real};
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
function enumerate_motif_embeddings(lattice::AbstractLattice, coords::AbstractVector{<:Tuple};
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

function _motif_geometry(lattice::MLattice)
    dims = lattice.supercell_dimensions
    scv = lattice.lattice_vectors * Diagonal([dims[1], dims[2], dims[3]])
    return scv, inv(scv), lattice.positions, lattice.periodicity
end

function _motif_geometry(lattice::AtomicLattice)
    _validate_atomic_components(lattice)
    _validate_atomic_ase_cache(lattice)
    cell = pyconvert(Matrix{Float64}, lattice.ase_lattice.get_cell())
    scv = permutedims(cell)
    singular_2d = !lattice.periodicity[3] && all(iszero, view(scv, :, 3))
    if singular_2d
        rec = zeros(3, 3)
        rec[1:2, 1:2] = inv(scv[1:2, 1:2])
    else
        rec = inv(scv)
    end
    return scv, rec, lattice.lattice_positions, lattice.periodicity
end

function _enumerate_motif_core(lattice::AbstractLattice, distances::Vector{Float64}, K::Int,
                               template::Union{Nothing,Matrix{Float64}},
                               tol::Float64, expected_count::Union{Int,Nothing})
    npairs = length(distances)
    sig = sort(distances)
    all(>(0.0), sig) || throw(ArgumentError("motif distances must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    M = num_sites(lattice)
    scv, rec, pos, per = _motif_geometry(lattice)

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

"""
    replicate_walkers(template::MLattice{C,G}, K::Int) -> Vector{LatticeWalker}

Build `K` walkers whose configurations share the template's run-invariant
geometry (lattice vectors, positions, basis, supercell dimensions,
periodicity, cutoff radii, neighbor lists, and adsorption mask) by
reference, each with its own independent occupancy vectors copied from the
template. The geometry fields are written only during construction and the
Monte Carlo kernels mutate occupancies exclusively, so sharing is safe on
the serial lattice drivers; relative to the `deepcopy(template)` idiom the
saving is one whole neighbor nest per walker. Walkers start at
`energy = 0.0u"eV"`, `iter = 0`, matching the deepcopy idiom's usual
construction.
"""
function replicate_walkers(template::MLattice{C,G}, K::Int) where {C,G}
    return [LatticeWalker(
                MLattice{C,G}(Val(:share_geometry), template,
                              [copy(v) for v in template.components]),
                energy=0.0u"eV", iter=0)
            for _ in 1:K]
end


"""Build `K` independent AtomicLattice walkers, including independent ASE caches."""
function replicate_walkers(template::AtomicLattice, K::Int)
    K >= 0 || throw(ArgumentError("K must be nonnegative, got $K"))
    return [LatticeWalker(deepcopy(template), energy=0.0u"eV", iter=0)
            for _ in 1:K]
end
