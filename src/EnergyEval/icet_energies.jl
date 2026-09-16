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

# ICET is not a dependency of FreeBird

`pyimport("icet")` happens inside the constructor below, never at module scope,
so `using FreeBird` neither imports nor requires ICET. Only constructing an
`ICETHamiltonian` does.

This is worth being explicit about, because the obvious alternative does not
work: a Julia package extension cannot help here. `[weakdeps]`/`[extensions]`
trigger on a **Julia** package being loaded, and ICET is a **Python** module
reached through PythonCall — there is nothing to weakly depend on. Keyed on
PythonCall, which is a hard dependency, an extension would load unconditionally
and achieve nothing. Deferring the import is the mechanism that actually gives
the property.
"""
struct ICETHamiltonian <: ClassicalHamiltonian
    calculator::Py
    n_sites::Int
    E_clean::Float64
    julia_to_icet::Vector{Int}
end

"""
    ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice)

Load a cluster expansion from `ce_path` and build a calculator for the supercell
of `base_lattice`.

Requires ICET and mchammer to be importable by the Python interpreter PythonCall
is using; the import happens here and nowhere else.

The energy of the empty lattice is measured once and stored as `E_clean`, so
`interacting_energy` returns adsorption energies referenced to the bare surface
rather than absolute cluster-expansion values.
"""
function ICETHamiltonian(ce_path::String, base_lattice::AtomicLattice)
    num_lattice_components(base_lattice) == 1 || throw(ArgumentError(
        "ICETHamiltonian represents a binary empty/occupied cluster " *
        "expansion and requires a one-species AtomicLattice"))
    base_lattice.periodicity[1:2] == (true, true) || throw(ArgumentError(
        "ICETHamiltonian repeats a periodic cluster-expansion model and requires " *
        "periodicity=(true, true, ...) on the AtomicLattice"))
    base_lattice.surface == :fcc100 || throw(ArgumentError(
        "ICETHamiltonian maps fcc(100) ICET models and " *
        "requires surface=:fcc100; got surface=$(base_lattice.surface)"))

    icet = pyimport("icet")
    mcham = pyimport("mchammer.calculators")

    ce = icet.ClusterExpansion.read(ce_path)
    ce_prim = ce._cluster_space.primitive_structure
    nx, ny = base_lattice.supercell_dimensions[1], base_lattice.supercell_dimensions[2]
    model = ce_prim.repeat(pytuple((nx, ny, 1)))
    calc = mcham.ClusterExpansionCalculator(model, ce)

    n_sites = nx * ny
    if length(base_lattice.all_sites) != n_sites
        throw(ArgumentError(
            "ICET model has $n_sites sites (from supercell $(nx)x$(ny)) but the " *
            "lattice has $(length(base_lattice.all_sites)) adsorption sites. " *
            "A cluster expansion over one site type cannot be mapped onto a " *
            "lattice built with type_of_sites = $(base_lattice.type_of_sites)."))
    end

    # ICET reports meV; every energy in FreeBird is eV.
    occ_empty = pylist(zeros(Int, n_sites))
    E_clean = pyconvert(Float64, calc.calculate_total(occupations=occ_empty)) / 1000

    icet_pos = pyconvert(Matrix{Float64}, model.get_positions())[:, 1:2]
    icet_to_julia = build_icet_to_julia_map(base_lattice, icet_pos)
    julia_to_icet = zeros(Int, n_sites)
    for (icet_idx, julia_idx) in enumerate(icet_to_julia)
        julia_to_icet[julia_idx] = icet_idx
    end

    @debug "ICETHamiltonian built" n_sites E_clean
    return ICETHamiltonian(calc, n_sites, E_clean, julia_to_icet)
end

"""
    build_icet_to_julia_map(base_lattice::AtomicLattice, icet_positions::Matrix{Float64})

Match ICET's site ordering to `base_lattice.all_sites`, returning a vector that
maps an ICET index to the corresponding index into `all_sites`.

The two orderings have no reason to agree, and the coordinate systems can differ
by a rigid translation, so the offset is measured from the first ICET site
before matching. Distances are compared under the minimum-image convention in
the surface plane.

Throws if any ICET site has no partner, and asserts the result is a bijection —
a silently non-injective map would send two ICET sites to one lattice site and
quietly lose occupancy.
"""
function build_icet_to_julia_map(base_lattice::AtomicLattice, icet_positions::Matrix{Float64})
    n_sites = length(base_lattice.all_sites)
    tol = nn_distance(base_lattice) * 0.1
    # `ase.build.fcc100` uses the surface nearest-neighbour spacing a/sqrt(2),
    # not the cubic lattice constant `a`, in the two in-plane cell vectors.
    # Read those vectors from the actual slab: using `a * nx/ny` here gives a
    # plausible but wrong minimum-image map near a periodic boundary.
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

    @assert length(unique(icet_to_julia)) == n_sites "ICET to Julia map is not bijective"
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
    occ = zeros(Int, h.n_sites)
    num_lattice_components(lattice) == 1 || throw(ArgumentError(
        "ICETHamiltonian requires a one-species AtomicLattice"))
    for julia_idx in occupied_indices(lattice)
        occ[h.julia_to_icet[julia_idx]] = 1
    end
    E_meV = pyconvert(Float64, h.calculator.calculate_total(occupations=pylist(occ)))
    return (E_meV / 1000 - h.E_clean) * u"eV"
end
