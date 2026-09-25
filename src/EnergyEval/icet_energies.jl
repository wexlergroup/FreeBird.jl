"""
    ICETHamiltonian <: ClassicalHamiltonian

A cluster-expansion Hamiltonian evaluated by ICET's `ClusterExpansionCalculator`.

# Fields
- `calculator::Py`: the `mchammer.calculators.ClusterExpansionCalculator`.
- `n_sites::Int`: number of sites the ICET model was built for.
- `E_clean::Float64`: energy of the empty lattice, in eV. Subtracted from every
  evaluation so that energies are referenced to the bare surface.
- `julia_to_icet::Vector{Int}`: site permutation, mapping an index into
  `all_sites` to the corresponding ICET site index.
- `empty_occupation::Int` and `occupied_occupation::Int`: ICET atomic-number
  codes for an empty and occupied adsorption site.
- `energy_scale_to_eV::Float64`: multiplier converting the cluster expansion's
  modeled-property unit to eV.
- `adsorbate_atom::String`: adsorbate species used when the Hamiltonian was
  constructed.
- `geometry_signature::NamedTuple`: surface, cell, boundary, site-family, and
  ordered-site identity of the lattice used to build the calculator mapping.

# ICET is not a dependency of FreeBird

`pyimport("icet")` happens inside the constructor below, never at module scope,
so `using FreeBird` neither imports nor requires ICET. Only constructing an
`ICETHamiltonian` does.

ICET is a Python module and cannot activate a Julia package extension.
Importing it inside the constructor keeps ICET optional when loading FreeBird.
"""
struct ICETHamiltonian <: ClassicalHamiltonian
    calculator::Py
    n_sites::Int
    E_clean::Float64
    julia_to_icet::Vector{Int}
    empty_occupation::Int
    occupied_occupation::Int
    energy_scale_to_eV::Float64
    adsorbate_atom::String
    geometry_signature::NamedTuple
end

_icet_geometry_signature(lattice::AtomicLattice) = (
    lattice_atom=lattice.lattice_atom,
    surface=lattice.surface,
    dimensions=lattice.supercell_dimensions,
    periodicity=lattice.periodicity,
    lattice_constant=lattice.lattice_constant,
    lattice_constant_c=lattice.lattice_constant_c,
    type_of_sites=Tuple(lattice.type_of_sites),
    all_sites=Tuple(lattice.all_sites),
)

"""
    ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice;
                    empty_occupation=nothing,
                    occupied_occupation=nothing,
                    energy_scale_to_eV=1.0)

Load a cluster expansion from `ce_path` and build a calculator for the supercell
of `base_lattice`.

Requires ICET and mchammer to be importable by the Python interpreter PythonCall
is using; the import happens here and nowhere else.

`empty_occupation` and `occupied_occupation` are the atomic-number codes used
by the binary ICET model. The occupied code defaults to the atomic number of
the lattice's sole adsorbate species; the empty code is inferred from ICET's
allowed occupations when that inference is unique. Pass either keyword
explicitly for models using different pseudo-species codes.

`energy_scale_to_eV` converts the modeled property returned by ICET to eV. Its
default, `1.0`, assumes that the cluster expansion already returns eV; use
`1e-3` for a model expressed in meV.

The energy of the empty lattice is measured once and stored as `E_clean`, so
`interacting_energy` returns adsorption energies referenced to the bare surface
rather than absolute cluster-expansion values.
"""
function ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice;
                         empty_occupation::Union{Nothing,Integer}=nothing,
                         occupied_occupation::Union{Nothing,Integer}=nothing,
                         energy_scale_to_eV::Real=1.0)
    AbstractWalkers._validate_atomic_components(base_lattice)
    num_lattice_components(base_lattice) == 1 || throw(ArgumentError(
        "ICETHamiltonian represents a binary empty/occupied cluster " *
        "expansion and requires a one-species AtomicLattice"))
    base_lattice.periodicity[1:2] == (true, true) || throw(ArgumentError(
        "ICETHamiltonian repeats a periodic cluster-expansion model and requires " *
        "periodicity=(true, true, ...) on the AtomicLattice"))
    base_lattice.surface == :fcc100 || throw(ArgumentError(
        "ICETHamiltonian maps fcc(100) ICET models and " *
        "requires surface=:fcc100; got surface=$(base_lattice.surface)"))
    isfinite(energy_scale_to_eV) && energy_scale_to_eV > 0 || throw(ArgumentError(
        "energy_scale_to_eV must be finite and positive, got $energy_scale_to_eV"))

    icet = pyimport("icet")
    mcham = pyimport("mchammer.calculators")
    ase_data = pyimport("ase.data")

    ce = icet.ClusterExpansion.read(ce_path)
    ce_prim = ce.get_cluster_space_copy().primitive_structure
    nx, ny = base_lattice.supercell_dimensions[1], base_lattice.supercell_dimensions[2]
    model = ce_prim.repeat(pytuple((nx, ny, 1)))
    _validate_icet_model_cell(
        base_lattice,
        pyconvert(Matrix{Float64}, model.get_cell()),
        Tuple(pyconvert(Vector{Bool}, model.get_pbc())))
    n_sites = length(model)
    if length(base_lattice.all_sites) != n_sites
        throw(ArgumentError(
            "the repeated ICET model has $n_sites sites but the " *
            "lattice has $(length(base_lattice.all_sites)) adsorption sites. " *
            "A cluster expansion over one site type cannot be mapped onto a " *
            "lattice built with type_of_sites = $(base_lattice.type_of_sites)."))
    end
    calc = mcham.ClusterExpansionCalculator(model, ce)

    occupied_code = occupied_occupation === nothing ? pyconvert(
        Int, ase_data.atomic_numbers[only(base_lattice.adsorbate_atoms)]) :
        Int(occupied_occupation)
    allowed_by_site = [sort!(pyconvert(Vector{Int},
        calc.sublattices.get_allowed_numbers_on_site(i))) for i in 0:(n_sites - 1)]
    empty_code = if empty_occupation === nothing
        alternatives = unique(vcat(
            [filter(!=(occupied_code), allowed) for allowed in allowed_by_site]...))
        length(alternatives) == 1 || throw(ArgumentError(
            "could not infer one empty-site atomic-number code from ICET's " *
            "allowed occupations $(allowed_by_site); pass empty_occupation explicitly"))
        only(alternatives)
    else
        Int(empty_occupation)
    end
    empty_code != occupied_code || throw(ArgumentError(
        "empty_occupation and occupied_occupation must differ"))
    expected_codes = sort([empty_code, occupied_code])
    all(==(expected_codes), allowed_by_site) || throw(ArgumentError(
        "ICETHamiltonian requires every ICET site to allow exactly the empty/occupied " *
        "atomic numbers $expected_codes; got $(allowed_by_site)"))

    occ_empty = pylist(fill(empty_code, n_sites))
    E_clean = pyconvert(Float64, calc.calculate_total(occupations=occ_empty)) *
              Float64(energy_scale_to_eV)

    icet_pos = pyconvert(Matrix{Float64}, model.get_positions())[:, 1:2]
    icet_to_julia = build_icet_to_julia_map(base_lattice, icet_pos)
    julia_to_icet = zeros(Int, n_sites)
    for (icet_idx, julia_idx) in enumerate(icet_to_julia)
        julia_to_icet[julia_idx] = icet_idx
    end

    @debug "ICETHamiltonian built" n_sites E_clean
    return ICETHamiltonian(calc, n_sites, E_clean, julia_to_icet,
                           empty_code, occupied_code,
                           Float64(energy_scale_to_eV),
                           only(base_lattice.adsorbate_atoms),
                           _icet_geometry_signature(base_lattice))
end

"""Validate that the repeated ICET model and adsorption lattice share a metric."""
function _validate_icet_model_cell(base_lattice::AtomicLattice,
                                   model_cell::AbstractMatrix{<:Real},
                                   model_pbc::Tuple{Bool,Bool,Bool})
    model_pbc == base_lattice.periodicity || throw(ArgumentError(
        "the repeated ICET model periodicity $model_pbc does not match the " *
        "AtomicLattice periodicity $(base_lattice.periodicity)"))
    size(model_cell, 1) >= 2 && size(model_cell, 2) >= 3 ||
        throw(DimensionMismatch(
            "the ICET model cell must contain two three-dimensional surface vectors"))
    AbstractWalkers._validate_atomic_ase_cache(base_lattice)
    lattice_cell = pyconvert(Matrix{Float64}, base_lattice.ase_lattice.get_cell())
    model_vectors = Matrix{Float64}(model_cell[1:2, 1:3])
    lattice_vectors = lattice_cell[1:2, 1:3]
    model_metric = model_vectors * transpose(model_vectors)
    lattice_metric = lattice_vectors * transpose(lattice_vectors)
    isapprox(model_metric, lattice_metric; rtol=1e-8, atol=1e-10) ||
        throw(ArgumentError(
            "the repeated ICET model cell does not match the AtomicLattice " *
            "in-plane metric; verify the lattice constant and supercell"))
    return nothing
end

"""
    build_icet_to_julia_map(base_lattice::AtomicLattice, icet_positions::Matrix{Float64})

Match ICET's site ordering to `base_lattice.all_sites`, returning a vector that
maps an ICET index to the corresponding index into `all_sites`.

The two orderings have no reason to agree, and the coordinate systems can differ
by a rigid translation, so the offset is measured from the first ICET site
before matching. Distances are compared under the minimum-image convention in
the surface plane.

Throws if any ICET site has no partner or if the result is not a bijection.
"""
function build_icet_to_julia_map(base_lattice::AtomicLattice, icet_positions::Matrix{Float64})
    n_sites = length(base_lattice.all_sites)
    tol = nn_distance(base_lattice) * 0.1
    # `ase.build.fcc100` uses the surface nearest-neighbour spacing a/sqrt(2),
    # not the cubic lattice constant `a`, in the two in-plane cell vectors.
    # Read those vectors from the actual slab: using `a * nx/ny` here gives a
    # plausible but wrong minimum-image map near a periodic boundary.
    AbstractWalkers._validate_atomic_ase_cache(base_lattice)
    cell = pyconvert(Matrix{Float64}, base_lattice.ase_lattice.get_cell())
    cell_x = hypot(cell[1, 1], cell[1, 2])
    cell_y = hypot(cell[2, 1], cell[2, 2])

    # Minimum-image separation in the surface plane.
    function wrapped_delta(a, b, cell)
        d = mod(abs(a - b), cell)
        return min(d, cell - d)
    end

    # The two coordinate systems can differ by a rigid translation. Measure it
    # from the first ICET site and its nearest lattice site.
    icet_x0, icet_y0 = icet_positions[1, 1], icet_positions[1, 2]
    best_dist = Inf
    offset_x = 0.0
    offset_y = 0.0
    for s in base_lattice.all_sites
        d = sqrt(wrapped_delta(s[1], icet_x0, cell_x)^2 +
                 wrapped_delta(s[2], icet_y0, cell_y)^2)
        if d < best_dist
            best_dist = d
            offset_x = s[1] - icet_x0
            offset_y = s[2] - icet_y0
        end
    end
    @debug "ICET to Julia coordinate offset" offset_x offset_y

    icet_to_julia = zeros(Int, n_sites)
    for icet_idx in 1:n_sites
        icet_x = icet_positions[icet_idx, 1] + offset_x
        icet_y = icet_positions[icet_idx, 2] + offset_y
        julia_idx = findfirst(base_lattice.all_sites) do s
            wrapped_delta(s[1], icet_x, cell_x) < tol &&
                wrapped_delta(s[2], icet_y, cell_y) < tol
        end
        julia_idx === nothing && error(
            "no lattice site within $tol A of ICET site $(icet_idx - 1) " *
            "at ($icet_x, $icet_y) after the offset")
        icet_to_julia[icet_idx] = julia_idx
    end

    length(unique(icet_to_julia)) == n_sites || throw(ArgumentError(
        "ICET to Julia site map is not bijective: $icet_to_julia"))
    return icet_to_julia
end

"""
    interacting_energy(lattice::AtomicLattice, h::ICETHamiltonian)

Cluster-expansion energy of the lattice's current occupancy, in eV, referenced
to the empty lattice.

The result is a `Unitful` energy. Occupied site indices map directly into ICET
ordering through `julia_to_icet`, making occupancy assembly linear in the
number of lattice sites.
"""
function interacting_energy(lattice::AtomicLattice, h::ICETHamiltonian)
    AbstractWalkers._validate_atomic_components(lattice)
    num_lattice_components(lattice) == 1 || throw(ArgumentError(
        "ICETHamiltonian requires a one-species AtomicLattice"))
    num_sites(lattice) == h.n_sites || throw(DimensionMismatch(
        "ICETHamiltonian was built for $(h.n_sites) sites, but the lattice has " *
        "$(num_sites(lattice))"))
    length(h.julia_to_icet) == h.n_sites || throw(DimensionMismatch(
        "ICETHamiltonian site map has $(length(h.julia_to_icet)) entries, " *
        "expected $(h.n_sites)"))
    only(lattice.adsorbate_atoms) == h.adsorbate_atom || throw(ArgumentError(
        "ICETHamiltonian was built for adsorbate $(h.adsorbate_atom), but the " *
        "lattice uses $(only(lattice.adsorbate_atoms))"))
    _icet_geometry_signature(lattice) == h.geometry_signature ||
        throw(ArgumentError(
            "ICETHamiltonian can only evaluate the AtomicLattice geometry used " *
            "to build its site mapping"))
    occ = fill(h.empty_occupation, h.n_sites)
    for julia_idx in occupied_indices(lattice)
        occ[h.julia_to_icet[julia_idx]] = h.occupied_occupation
    end
    raw_energy = pyconvert(Float64,
        h.calculator.calculate_total(occupations=pylist(occ)))
    return (raw_energy * h.energy_scale_to_eV - h.E_clean) * u"eV"
end
