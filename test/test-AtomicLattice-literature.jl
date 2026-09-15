using LinearAlgebra
using Test
using FreeBird

# These tests intentionally validate quantities that AtomicLattice can
# reproduce without a fitted potential: crystallographic site geometry,
# overlayer coverage/periodicity, and measured bond geometry.  Adsorption
# energies and transition temperatures belong in a separate model-validation
# suite because they require the paper's Hamiltonian or an equivalent MLIP.

const _LIT_AW = FreeBird.AbstractWalkers

"""Integer primitive-cell coordinates, independent of the site-list order."""
function _lit_grid_indices(lattice::AtomicLattice)
    d1, d2, _ = lattice.supercell_dimensions
    fractions = _LIT_AW._atomic_site_fractions(lattice)
    f1zero = minimum(f[1] for f in fractions)
    f2zero = minimum(f[2] for f in fractions)
    return [(mod(round(Int, d1 * (f[1] - f1zero)), d1),
             mod(round(Int, d2 * (f[2] - f2zero)), d2))
            for f in fractions]
end

function _lit_configuration(base::AtomicLattice, predicate)
    lattice = deepcopy(base)
    indices = _lit_grid_indices(lattice)
    lattice.components[1] .= [predicate(i, j) for (i, j) in indices]
    lattice.ase_dirty = true
    return lattice
end

function _lit_cell_geometry(lattice::AtomicLattice)
    cell = _LIT_AW.pyconvert(Matrix{Float64}, lattice.ase_lattice.get_cell())
    inplane = [cell[1, 1] cell[2, 1]; cell[1, 2] cell[2, 2]]
    return inplane, inv(inplane)
end

function _lit_minimum_occupied_distance(lattice::AtomicLattice)
    occupied = occupied_indices(lattice)
    fractions = _LIT_AW._atomic_site_fractions(lattice)
    inplane, _ = _lit_cell_geometry(lattice)
    distances = Float64[]
    for a in 1:(length(occupied) - 1), b in (a + 1):length(occupied)
        delta = fractions[occupied[b]] - fractions[occupied[a]]
        delta .-= round.(delta)
        push!(distances, norm(inplane * delta))
    end
    return minimum(distances)
end

@testset "AtomicLattice literature validation" begin
    @testset "O/Pd(100) ordered hollow-site overlayers" begin
        # Zhang, Blum & Reuter, Phys. Rev. B 75, 235406 (2007),
        # doi:10.1103/PhysRevB.75.235406.  The paper uses a=3.947 Angstrom,
        # puts O on fourfold Pd(100) hollow sites, and reports the ordered
        # coverages 1/9, 1/4 [p(2x2)], 1/2 [c(2x2)], 3/4, and 1 ML.
        base = AtomicLattice{1,SquareLattice}(
            lattice_atom="Pd", surface=:fcc100,
            supercell_dimensions=(6, 6, 5), lattice_constant=3.947,
            periodicity=(true, true, false), adsorbate_atoms=["O"],
            components=[fill(false, 36)], num_nearest_neighbors=2,
            type_of_sites=["hollow"])

        one_ninth = _lit_configuration(base, (i, j) -> i % 3 == 0 && j % 3 == 0)
        p2x2 = _lit_configuration(base, (i, j) -> iseven(i) && iseven(j))
        c2x2 = _lit_configuration(base, (i, j) -> iseven(i + j))
        three_quarters = _lit_configuration(base, (i, j) -> !(iseven(i) && iseven(j)))
        full = _lit_configuration(base, (i, j) -> true)

        @test collect(coverage.((one_ninth, p2x2, c2x2, three_quarters, full))) ≈
              [1 / 9, 1 / 4, 1 / 2, 3 / 4, 1.0]
        @test order_parameter_c2x2(c2x2) ≈ 1 / 2 atol=1e-12
        @test bragg_amplitude(p2x2, 3, 0) ≈ 1 / 4 atol=1e-12
        @test bragg_amplitude(p2x2, 0, 3) ≈ 1 / 4 atol=1e-12
        @test bragg_amplitude(p2x2, 3, 3) ≈ 1 / 4 atol=1e-12

        # The paper gives 8.37 Angstrom as the shortest O-O separation in the
        # sparse (3x3)-O structure.  This follows directly from the generated
        # adsorption lattice, rather than being inserted as a model parameter.
        @test _lit_minimum_occupied_distance(one_ninth) ≈ 8.37 atol=0.01

        # One hollow-site orbit per surface Pd atom and four equal nearest
        # top-layer neighbours verify the lattice and "fourfold hollow" parts
        # of the experimental model.
        inplane, reciprocal = _lit_cell_geometry(base)
        area = abs(det(inplane))
        @test num_sites(base) / area * 1e16 ≈ 2 / 3.947^2 * 1e16 rtol=1e-12

        slab_positions = _LIT_AW.pyconvert(
            Matrix{Float64}, base.ase_lattice.get_positions())
        top_z = maximum(slab_positions[:, 3])
        top = [view(slab_positions, i, 1:2) for i in axes(slab_positions, 1)
               if isapprox(slab_positions[i, 3], top_z; atol=1e-10)]
        hollow_fraction = reciprocal * collect(first(base.all_sites))
        top_distances = sort(map(top) do position
            delta = reciprocal * collect(position) - hollow_fraction
            delta .-= round.(delta)
            norm(inplane * delta)
        end)
        @test top_distances[1:4] ≈ fill(3.947 / 2, 4) atol=1e-10
        @test top_distances[5] > top_distances[4]
    end

    @testset "O/Pd(111) p(2x2) overlayer" begin
        # Zheng & Altman, Surf. Sci. 462, 151-168 (2000),
        # doi:10.1016/S0039-6028(00)00599-9, report a (2x2) O structure at
        # 0.25 ML on Pd(111).  A single fcc-site orbit must represent that
        # coverage and its three symmetry-equivalent M-point peaks.
        base = AtomicLattice{1,TriangularLattice}(
            lattice_atom="Pd", surface=:fcc111,
            supercell_dimensions=(4, 4, 4), lattice_constant=3.947,
            periodicity=(true, true, false), adsorbate_atoms=["O"],
            components=[fill(false, 16)], num_nearest_neighbors=1,
            type_of_sites=["fcc"])
        p2x2 = _lit_configuration(base, (i, j) -> iseven(i) && iseven(j))

        @test coverage(p2x2) == 1 / 4
        @test order_parameter_p2x2(p2x2) ≈ sqrt(3) / 4 atol=1e-12
        @test bragg_amplitude(p2x2, 2, 0) ≈ 1 / 4 atol=1e-12
        @test bragg_amplitude(p2x2, 0, 2) ≈ 1 / 4 atol=1e-12
        @test bragg_amplitude(p2x2, 2, 2) ≈ 1 / 4 atol=1e-12
    end

    @testset "Cl/Ni(111) sqrt(3) adsorption geometry" begin
        # Wang et al., Phys. Rev. B 44, 13711-13719 (1991),
        # doi:10.1103/PhysRevB.44.13711, find Cl in the fcc threefold hollow
        # site, 1.837(8) Angstrom above the first Ni layer, with a 2.332(6)
        # Angstrom Cl-Ni bond.  a=3.5238 Angstrom is the room-temperature fcc
        # Ni lattice constant tabulated in NBS Circular 592.
        base = AtomicLattice{1,TriangularLattice}(
            lattice_atom="Ni", surface=:fcc111,
            supercell_dimensions=(3, 3, 4), lattice_constant=3.5238,
            periodicity=(true, true, false), adsorbate_atoms=["Cl"],
            components=[fill(false, 9)], num_nearest_neighbors=1,
            type_of_sites=["fcc"], adsorbate_height=1.837)
        sqrt3 = _lit_configuration(base, (i, j) -> mod(i - j, 3) == 0)
        sync_ase_lattice!(sqrt3)

        @test coverage(sqrt3) == 1 / 3
        @test order_parameter_sqrt3(sqrt3) ≈ 1 / 3 atol=1e-12

        slab = sqrt3.ase_lattice
        tags = _LIT_AW.pyconvert(Vector{Int}, slab.get_tags())
        positions = _LIT_AW.pyconvert(Matrix{Float64}, slab.get_positions())
        adsorbates = findall(==(0), tags)
        substrate = findall(!=(0), tags)
        top_z = maximum(positions[substrate, 3])
        @test all(i -> isapprox(positions[i, 3] - top_z, 1.837; atol=1e-12),
                  adsorbates)

        for adsorbate in adsorbates
            distances = sort(_LIT_AW.pyconvert(Vector{Float64},
                slab.get_distances(adsorbate - 1, substrate .- 1; mic=true)))
            @test length(filter(d -> abs(d - 2.332) <= 0.006, distances)) == 3
            @test distances[4] > 2.332 + 0.006
        end
    end
end
