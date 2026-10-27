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

ICET is optional and imported only when constructing an `ICETHamiltonian`.
Loading FreeBird does not require ICET.
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

function _resolve_icet_occupation_codes(
        allowed_by_site::AbstractVector{<:AbstractVector{<:Integer}},
        adsorbate_code::Integer,
        empty_occupation::Union{Nothing,Integer},
        occupied_occupation::Union{Nothing,Integer})
    isempty(allowed_by_site) && throw(ArgumentError(
        "the ICET model contains no occupation sites"))
    allowed = unique([sort!(Int.(collect(codes))) for codes in allowed_by_site])
    length(allowed) == 1 || throw(ArgumentError(
        "ICETHamiltonian requires the same binary occupations on every site; " *
        "got $(allowed_by_site)"))
    codes = only(allowed)
    length(codes) == 2 || throw(ArgumentError(
        "ICETHamiltonian requires exactly two occupation codes per site; " *
        "got $codes"))

    empty_code = empty_occupation === nothing ? nothing : Int(empty_occupation)
    occupied_code = occupied_occupation === nothing ? nothing : Int(occupied_occupation)

    if occupied_code === nothing
        if Int(adsorbate_code) in codes
            occupied_code = Int(adsorbate_code)
        elseif empty_code !== nothing
            alternatives = filter(!=(empty_code), codes)
            length(alternatives) == 1 || throw(ArgumentError(
                "could not infer occupied_occupation from ICET occupation codes " *
                "$codes and empty_occupation=$empty_code"))
            occupied_code = only(alternatives)
        elseif 0 in codes
            # ASE/ICET use atomic number zero (`X`) for a vacancy. This also
            # supports models that use a pseudo-species for the adsorbate.
            occupied_code = only(filter(!=(0), codes))
        else
            throw(ArgumentError(
                "the adsorbate atomic-number code $adsorbate_code is not present " *
                "in ICET occupation codes $codes; pass occupied_occupation and " *
                "empty_occupation explicitly"))
        end
    end

    if empty_code === nothing
        alternatives = filter(!=(occupied_code), codes)
        length(alternatives) == 1 || throw(ArgumentError(
            "could not infer empty_occupation from ICET occupation codes $codes " *
            "and occupied_occupation=$occupied_code"))
        empty_code = only(alternatives)
    end

    empty_code != occupied_code || throw(ArgumentError(
        "empty_occupation and occupied_occupation must differ"))
    sort([empty_code, occupied_code]) == codes || throw(ArgumentError(
        "ICETHamiltonian occupation codes $(sort([empty_code, occupied_code])) " *
        "do not match the model's allowed codes $codes"))
    return empty_code, occupied_code
end

"""
    ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice;
                    empty_occupation=nothing,
                    occupied_occupation=nothing,
                    energy_scale_to_eV)

Load a cluster expansion from `ce_path` and build a calculator for the supercell
of `base_lattice`.

Requires a one-species, in-plane-periodic fcc(100) adsorption lattice. The ICET
primitive structure is repeated by `(nx, ny, 1)`; its in-plane cell and site
count must match the adsorption lattice. Incompatible models are rejected.

Requires ICET and mchammer to be importable by the Python interpreter PythonCall
is using; the import happens here and nowhere else.

`empty_occupation` and `occupied_occupation` are the atomic-number codes used
by the binary ICET model. The occupied code defaults to the atomic number of
the lattice's sole adsorbate species. A binary model using ASE's vacancy code
`X` (atomic number zero) and one pseudo-species is inferred automatically.
Pass both keywords explicitly for other pseudo-species conventions.

`energy_scale_to_eV` converts the modeled property returned by ICET to eV. Its
value must be supplied because ICET cluster-expansion files do not define a
standard energy unit. Use `1.0` for a model expressed in eV and `1e-3` for a
model expressed in meV.

The energy of the empty lattice is measured once and stored as `E_clean`, so
`interacting_energy` returns adsorption energies referenced to the bare surface
rather than absolute cluster-expansion values.
"""
function ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice;
                         empty_occupation::Union{Nothing,Integer}=nothing,
                         occupied_occupation::Union{Nothing,Integer}=nothing,
                         energy_scale_to_eV::Union{Nothing,Real}=nothing)
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
    energy_scale_to_eV === nothing && throw(ArgumentError(
        "energy_scale_to_eV must be supplied explicitly because ICET files do " *
        "not define a standard energy unit; use 1.0 for eV or 1e-3 for meV"))
    energy_scale = Float64(energy_scale_to_eV)
    isfinite(energy_scale) && energy_scale > 0 || throw(ArgumentError(
        "energy_scale_to_eV must remain finite and positive after conversion " *
        "to Float64, got $energy_scale_to_eV"))

    icet = pyimport("icet")
    mcham = pyimport("mchammer.calculators")
    ase_data = pyimport("ase.data")

    ce = icet.ClusterExpansion.read(ce_path)
    ce_prim = ce.get_cluster_space_copy().primitive_structure
    nx, ny = base_lattice.supercell_dimensions[1], base_lattice.supercell_dimensions[2]
    model = ce_prim.repeat(pytuple((nx, ny, 1)))
    source_pbc = Tuple(pyconvert(Vector{Bool}, model.get_pbc()))
    source_pbc[1:2] == base_lattice.periodicity[1:2] || throw(ArgumentError(
        "the repeated ICET model in-plane periodicity $(source_pbc[1:2]) does " *
        "not match the AtomicLattice periodicity $(base_lattice.periodicity[1:2])"))
    # ICET cluster spaces are three-dimensionally periodic, while an adsorption
    # lattice has a nonperiodic surface-normal direction. The cluster basis and
    # site mapping use only the in-plane cell, so adapt z before constructing
    # the calculator and retain strict validation of x and y above.
    model.set_pbc(pytuple(base_lattice.periodicity))
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

    allowed_by_site = [sort!(pyconvert(Vector{Int},
        calc.sublattices.get_allowed_numbers_on_site(i))) for i in 0:(n_sites - 1)]
    adsorbate_code = pyconvert(
        Int, ase_data.atomic_numbers[only(base_lattice.adsorbate_atoms)])
    empty_code, occupied_code = _resolve_icet_occupation_codes(
        allowed_by_site, adsorbate_code, empty_occupation, occupied_occupation)

    occ_empty = pylist(fill(empty_code, n_sites))
    E_clean = pyconvert(Float64, calc.calculate_total(occupations=occ_empty)) *
              energy_scale
    isfinite(E_clean) || throw(DomainError(E_clean,
        "the ICET model must have a finite empty-lattice energy in eV"))

    icet_pos = pyconvert(Matrix{Float64}, model.get_positions())[:, 1:2]
    icet_to_julia = build_icet_to_julia_map(base_lattice, icet_pos)
    julia_to_icet = zeros(Int, n_sites)
    for (icet_idx, julia_idx) in enumerate(icet_to_julia)
        julia_to_icet[julia_idx] = icet_idx
    end

    @debug "ICETHamiltonian built" n_sites E_clean
    return ICETHamiltonian(calc, n_sites, E_clean, julia_to_icet,
                           empty_code, occupied_code,
                           energy_scale,
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
    sort(h.julia_to_icet) == collect(1:h.n_sites) || throw(ArgumentError(
        "ICETHamiltonian site map must be a permutation of 1:$(h.n_sites), " *
        "got $(h.julia_to_icet)"))
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
    energy = raw_energy * h.energy_scale_to_eV - h.E_clean
    isfinite(energy) || throw(DomainError(energy,
        "the ICET model must return a finite adsorption energy in eV"))
    return energy * u"eV"
end
