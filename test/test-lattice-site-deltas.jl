# Exactness coverage for the O(z) single-site flip deltas and the
# supports_site_deltas trait. All assertions are exact or fixed-tolerance;
# there are no statistical gates. The atol 1e-13 eV carries a 12x margin
# over the measured worst case of the prototype (<= 8e-15 eV over 200
# random flips per size).

# Minimal trait fixture: testset bodies cannot define structs, so it sits
# at the file's top level (identical re-definition on a double include is
# safe).
struct SFDMinimalHam <: FreeBird.AbstractHamiltonians.ClassicalHamiltonian end

# Trait-less evaluable fixture: a pair Hamiltonian wrapper that evaluates
# exactly as its base but does not declare `supports_site_deltas`, so an
# incremental walk under it must fall back to the full recompute. Every
# shipped single-component lattice Hamiltonian type now carries the trait,
# so the fallback path is reachable only through a fixture like this one.
struct SFDNoDeltaHam{H} <: FreeBird.AbstractHamiltonians.ClassicalHamiltonian
    base::H
end
FreeBird.EnergyEval.interacting_energy(lat::SLattice, h::SFDNoDeltaHam) =
    interacting_energy(lat, h.base)

@testset "Lattice site-flip deltas" begin
    using Random
    using Unitful

    function sfd_maxdev(lat, h, seed; nflips=200)
        Random.seed!(seed)
        maxdev = 0.0
        M = length(lat.components[1])
        for _ in 1:nflips
            s = rand(1:M)
            e0 = interacting_energy(lat, h)
            d = site_flip_delta(lat, h, s)
            lat.components[1][s] = !lat.components[1][s]
            e1 = interacting_energy(lat, h)
            lat.components[1][s] = !lat.components[1][s]
            maxdev = max(maxdev, abs(ustrip(u"eV", (e1 - e0) - d)))
        end
        return maxdev
    end

    function sfd_square(L; cutoffs=[1.1], adsorptions=:full,
                        image_multiplicity=false)
        MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(L, L, 1),
            periodicity=(true, true, false), cutoff_radii=cutoffs,
            components=[[false for _ in 1:L*L]], adsorptions=adsorptions,
            image_multiplicity=image_multiplicity)
    end

    function sfd_fill!(lat, seed, theta)
        Random.seed!(seed)
        for i in eachindex(lat.components[1])
            lat.components[1][i] = rand() < theta
        end
        return lat
    end

    h1 = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
    h2 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

    # Swap counterpart of sfd_maxdev: two sequentially composed flips, the
    # second delta on the intermediate configuration, against the full
    # energy difference of the swap (non-null pairs only)
    function sfd_swap_maxdev(lat, h, seed; nswaps=200)
        Random.seed!(seed)
        maxdev = 0.0
        occ = lat.components[1]
        M = length(occ)
        done = 0
        while done < nswaps
            i = rand(1:M)
            j = rand(1:M)
            occ[i] == occ[j] && continue
            e0 = interacting_energy(lat, h)
            d = site_flip_delta(lat, h, i)
            occ[i] = !occ[i]
            d += site_flip_delta(lat, h, j)
            occ[j] = !occ[j]
            e1 = interacting_energy(lat, h)
            occ[j] = !occ[j]
            occ[i] = !occ[i]
            maxdev = max(maxdev, abs(ustrip(u"eV", (e1 - e0) - d)))
            done += 1
        end
        return maxdev
    end

    # Two-site triangular cell, one layer: three trio figures (face,
    # linear, obtuse) over two pair shells. The linear trio's K d_max equals
    # the a1 circumference, so its enumeration warns and follows the torus
    # convention.
    function sfd_tri(nx, ny)
        MLattice{1,TriangularLattice}(supercell_dimensions=(nx, ny, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1, 1.8],
            components=[[false for _ in 1:2*nx*ny]], adsorptions=:full)
    end
    tri_cell = sfd_tri(6, 4)
    tri_face = enumerate_motif_embeddings(tri_cell, [1.0, 1.0, 1.0];
                                          expected_count=96)
    tri_lin = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
        tri_cell, [1.0, 1.0, 2.0]; expected_count=144)
    tri_obt = enumerate_motif_embeddings(tri_cell, [1.0, 1.0, sqrt(3)];
                                         expected_count=288)
    hc_tri = ClusterLatticeHamiltonian(h2,
        [ClusterInteraction(0.012u"eV", tri_face),
         ClusterInteraction(-0.007u"eV", tri_lin),
         ClusterInteraction(0.004u"eV", tri_obt)])

    # Offset-stacked three-layer triangular cell (layer k shifted by
    # k (1/2, sqrt(3)/6) in plane, spacing h, periodic third axis): shell 1
    # is the six in-plane neighbors, shell 2 the six adjacent-layer
    # neighbors at sqrt(1/3 + h^2). Two figures spanning two layers: a bond
    # with the adjacent-layer site above it, and a second-neighbor pair
    # with the adjacent-layer site above their midpoint. Three layers are
    # not a faithful quotient for either, so both enumerations warn.
    function sfd_offset(nx, ny, nz; h=1.0, image_multiplicity=true)
        M = 2 * nx * ny * nz
        MLattice{1,TriangularLattice}(
            [1.0 0.0 0.5; 0.0 sqrt(3) sqrt(3)/6; 0.0 0.0 h],
            [(0.0, 0.0, 0.0), (0.5, sqrt(3)/2, 0.0)],
            (nx, ny, nz), (true, true, true), [1.05, 1.2],
            [zeros(Bool, M)], ones(Bool, M);
            image_multiplicity=image_multiplicity)
    end
    hl = GenericLatticeHamiltonian(-0.04, [-0.01, -0.006], u"eV")
    off_cell = sfd_offset(4, 4, 3)
    off_bond = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
        off_cell, [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, sqrt(3)/6, 1.0)];
        expected_count=576)
    off_second = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
        off_cell, [(0.0, 0.0, 0.0), (0.0, sqrt(3), 0.0), (0.0, 2 * sqrt(3) / 3, 1.0)])
    hc_off = ClusterLatticeHamiltonian(hl,
        [ClusterInteraction(0.012u"eV", off_bond),
         ClusterInteraction(-0.007u"eV", off_second)])

    @testset "exactness against full-energy differences" begin
        # square, one shell, full adsorption
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_001, 0.5), h1,
                         93_001) <= 1e-13
        # square, two shells
        @test sfd_maxdev(sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]),
                                   92_002, 0.5), h2, 93_002) <= 1e-13
        # partial adsorption mask
        @test sfd_maxdev(sfd_fill!(sfd_square(4;
                                              adsorptions=[1, 2, 3, 5, 8, 13]),
                                   92_003, 0.3), h1, 93_003) <= 1e-13
        # triangular geometry
        tri = SLattice{TriangularLattice}(supercell_dimensions=(6, 4, 1),
            cutoff_radii=[1.1], components=[[1, 5, 10, 20, 30, 40]],
            adsorptions=:full)
        @test sfd_maxdev(sfd_fill!(tri, 92_004, 0.5), h1, 93_004) <= 1e-13
        # image-multiplicity small cell: self-image entries and duplicated
        # image bonds sit under the delta's j == site convention
        imcell = sfd_square(2; image_multiplicity=true)
        @test sfd_maxdev(sfd_fill!(imcell, 92_005, 0.5), h1, 93_005) <= 1e-13
        # site-field wrapper
        fld = collect(0.001 .* (1:16)) .* u"eV"
        hsf = SiteFieldLatticeHamiltonian(h1, fld)
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_006, 0.5), hsf,
                         93_006) <= 1e-13
        # multi-component Hamiltonian on a single-component lattice
        mlham = MLatticeHamiltonian(1, [h1])
        @test sfd_maxdev(sfd_fill!(sfd_square(4), 92_007, 0.5), mlham,
                         93_007) <= 1e-13
    end

    @testset "sign symmetry is exact" begin
        # The accumulator visits the same terms with the opposite sign, so
        # the back-flip delta is the exact floating-point negation
        lat = sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]), 92_008, 0.5)
        Random.seed!(93_008)
        for _ in 1:50
            s = rand(1:36)
            d_fwd = site_flip_delta(lat, h2, s)
            lat.components[1][s] = !lat.components[1][s]
            d_back = site_flip_delta(lat, h2, s)
            lat.components[1][s] = !lat.components[1][s]
            @test d_back == -d_fwd
        end
    end

    @testset "swap composition" begin
        # Two sequentially composed flips, the second evaluated on the
        # intermediate configuration, reproduce the full-energy difference
        # of the swap at the single-flip rounding class, including adjacent
        # origin-destination pairs
        lat = sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]), 92_009, 0.5)
        Random.seed!(93_009)
        maxdev = 0.0
        for trial in 1:100
            i = rand(1:36)
            # force adjacency on odd trials: destination = a first-shell
            # neighbor entry of the origin
            j = isodd(trial) ? lat.neighbors[i][1][1] : rand(1:36)
            e0 = interacting_energy(lat, h2)
            d = site_flip_delta(lat, h2, i)
            lat.components[1][i] = !lat.components[1][i]
            d += site_flip_delta(lat, h2, j)
            lat.components[1][j] = !lat.components[1][j]
            e1 = interacting_energy(lat, h2)
            # restore
            lat.components[1][j] = !lat.components[1][j]
            lat.components[1][i] = !lat.components[1][i]
            maxdev = max(maxdev, abs(ustrip(u"eV", (e1 - e0) - d)))
        end
        @test maxdev <= 1e-13
    end

    @testset "trait contract and shell guard" begin
        fld = collect(0.001 .* (1:16)) .* u"eV"
        clham = ClusterLatticeHamiltonian(h1,
            [ClusterInteraction(0.1u"eV", [(1, 2, 3)])])
        @test supports_site_deltas(h1)
        @test supports_site_deltas(MLatticeHamiltonian(1, [h1]))
        @test supports_site_deltas(SiteFieldLatticeHamiltonian(h1, fld))
        @test !supports_site_deltas(SFDMinimalHam())
        @test !supports_site_deltas(SFDNoDeltaHam(h1))
        # The cluster Hamiltonian opted in when ClusterInteraction gained
        # its incidence lists (these two asserts were the negative case)
        @test supports_site_deltas(clham)
        # the wrapper delegates to its base
        @test supports_site_deltas(SiteFieldLatticeHamiltonian(clham, fld))
        # the delta enforces the same shell-count rule as the full sweep
        lat1 = sfd_square(4)
        @test_throws ArgumentError site_flip_delta(lat1, h2, 1)
    end

    @testset "cluster and multilayer deltas" begin
        # Cluster Hamiltonian: 200 random single flips and 200 random swaps
        # against interacting_energy differences, on the single-layer
        # two-site cell with three trio figures and on the offset
        # three-layer cell with two-layer figures (a figure spanning layers
        # is counted from every one of its sites)
        @test sfd_maxdev(sfd_fill!(tri_cell, 92_010, 0.5), hc_tri,
                         93_010) <= 1e-12
        @test sfd_swap_maxdev(tri_cell, hc_tri, 93_011) <= 1e-12
        @test sfd_maxdev(sfd_fill!(off_cell, 92_012, 0.5), hc_off,
                         93_012) <= 1e-12
        @test sfd_swap_maxdev(off_cell, hc_off, 93_013) <= 1e-12

        # Pair Hamiltonian on the multilayer cell under both neighbor
        # conventions
        @test sfd_maxdev(off_cell, hl, 93_014) <= 1e-12
        @test sfd_swap_maxdev(off_cell, hl, 93_015) <= 1e-12
        off_min = sfd_fill!(sfd_offset(4, 4, 3; image_multiplicity=false),
                            92_016, 0.5)
        @test sfd_maxdev(off_min, hl, 93_016) <= 1e-12
        @test sfd_swap_maxdev(off_min, hl, 93_017) <= 1e-12

        # A one-cell-wide multilayer cell under image_multiplicity = true:
        # the a1 circumference (1.0) lies inside shell 1, so every site
        # carries two self-image entries at half weight per image
        off_self = sfd_fill!(sfd_offset(1, 2, 3), 92_018, 0.5)
        @test all(count(==(s), off_self.neighbors[s][1]) == 2
                  for s in eachindex(off_self.components[1]))
        @test sfd_maxdev(off_self, hl, 93_018) <= 1e-12
        @test sfd_swap_maxdev(off_self, hl, 93_019) <= 1e-12
    end

    @testset "incremental fixed-N walk" begin
        # Same-seed kernel A/B: a call that never mentions `incremental`
        # against `incremental = false`, digit for digit
        function sfd_walk(lat, h, seed, emax, inc; n=500, perturb=1e-9)
            Random.seed!(seed)
            l = deepcopy(lat)
            w = LatticeWalker(l, energy=interacting_energy(l, h), iter=0)
            _, r, _ = inc === nothing ?
                MC_random_walk!(n, w, h, emax; energy_perturb=perturb) :
                MC_random_walk!(n, w, h, emax; energy_perturb=perturb,
                                incremental=inc)
            return w.energy.val, r, copy(w.configuration.components[1])
        end
        sq = sfd_fill!(sfd_square(6; cutoffs=[1.1, 1.5]), 92_020, 0.5)
        e_sq = ustrip(u"eV", interacting_energy(sq, h2))
        for (lat, h, emax) in ((sq, h2, e_sq + 0.01),
                               (tri_cell, hc_tri,
                                ustrip(u"eV", interacting_energy(tri_cell, hc_tri)) + 0.01))
            a = sfd_walk(lat, h, 93_020, emax, nothing)
            b = sfd_walk(lat, h, 93_020, emax, false)
            @test a[1] == b[1]
            @test a[2] == b[2]
            @test a[3] == b[3]
        end

        # Anchor-drift weld: after a 10^4-step seeded incremental walk with
        # a zero perturbation the stored walker energy IS the raw anchor,
        # so it must agree with a from-scratch recompute to 1e-12 eV. The
        # ceiling sits above every energy so every proposal is accepted and
        # the anchor advances 10^4 times.
        function sfd_weld(lat, h, seed; n=10_000, emax=10.0, perturb=0.0)
            Random.seed!(seed)
            l = deepcopy(lat)
            w = LatticeWalker(l, energy=interacting_energy(l, h), iter=0)
            _, r, _ = MC_random_walk!(n, w, h, emax; energy_perturb=perturb,
                                      incremental=true)
            drift = abs(ustrip(u"eV", w.energy -
                               interacting_energy(w.configuration, h)))
            return drift, r
        end
        fld_off = collect(0.001 .* (1:96)) .* u"eV"
        @test sfd_weld(off_cell, hl, 93_021)[1] <= 1e-12
        @test sfd_weld(off_cell, hc_off, 93_022)[1] <= 1e-12
        @test sfd_weld(off_cell, SiteFieldLatticeHamiltonian(hc_off, fld_off),
                       93_023)[1] <= 1e-12
        # With a perturbation the stored energy carries the last accepted
        # perturbation, bounded by half its width
        @test sfd_weld(off_cell, hc_off, 93_024; perturb=1e-9)[1] <=
              0.5e-9 + 1e-12
        # A ceiling near the start energy exercises the revert path inside
        # the same weld
        e_off = ustrip(u"eV", interacting_energy(off_cell, hc_off))
        drift_cold, rate_cold = sfd_weld(off_cell, hc_off, 93_025;
                                         emax=e_off + 0.02)
        @test drift_cold <= 1e-12
        @test 0.0 < rate_cold < 1.0
        # Null-pair dominance: a nearly full lattice proposes mostly
        # equal-occupancy pairs, which must contribute exactly zero
        full_cell = sfd_offset(4, 4, 3)
        full_cell.components[1] .= true
        full_cell.components[1][7] = false
        @test sfd_weld(full_cell, hc_off, 93_026; n=2000)[1] <= 1e-12

        # Trait fallback: an incremental = true walk under a Hamiltonian
        # without the trait matches the same-seed default digit for digit,
        # in the fixed-N kernel and in the grand-canonical kernel (whose
        # shipped fallback test used the cluster Hamiltonian before it
        # opted in)
        hnd = SFDNoDeltaHam(h2)
        a = sfd_walk(sq, hnd, 93_027, e_sq + 0.01, true)
        b = sfd_walk(sq, hnd, 93_027, e_sq + 0.01, false)
        @test a[1] == b[1]
        @test a[2] == b[2]
        @test a[3] == b[3]
        function sfd_gc_walk(inc)
            Random.seed!(93_028)
            l = deepcopy(sq)
            w = LatticeWalker(l, energy=interacting_energy(l, hnd), iter=0)
            _, r, _, _, _, _ = MC_grand_canonical_walk!(500, w, hnd, 1.0e3,
                0.0; p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
                incremental=inc)
            return w.energy.val, sum(w.configuration.components[1]), r
        end
        @test sfd_gc_walk(true) == sfd_gc_walk(false)

        # Multi-component walkers keep the full recompute: same-seed
        # identity under incremental = true
        ml = MLattice{2,SquareLattice}(supercell_dimensions=(4, 4, 1),
            components=[[1, 3, 6, 8, 11], [2, 4, 9, 13]])
        mlh = MLatticeHamiltonian(2, [h1, GenericLatticeHamiltonian(-0.02, [-0.005], u"eV"),
                                      GenericLatticeHamiltonian(-0.03, [-0.008], u"eV")])
        function sfd_ml_walk(inc)
            Random.seed!(93_029)
            l = deepcopy(ml)
            w = LatticeWalker(l, energy=interacting_energy(l, mlh), iter=0)
            _, r, _ = MC_random_walk!(300, w, mlh, 0.0; energy_perturb=1e-9,
                                      incremental=inc)
            return w.energy.val, r, deepcopy(w.configuration.components)
        end
        @test sfd_ml_walk(true) == sfd_ml_walk(false)
    end
end
