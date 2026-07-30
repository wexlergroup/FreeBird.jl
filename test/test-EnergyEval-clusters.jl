@testset "Cluster lattice Hamiltonian tests" begin
    using Random

    cluster_square(L; cutoffs=[1.1, 1.5]) = MLattice{1,SquareLattice}(
        lattice_constant=1.0,
        basis=[(0.0, 0.0, 0.0)],
        supercell_dimensions=(L, L, 1),
        periodicity=(true, true, false),
        cutoff_radii=cutoffs,
        components=[[false for _ in 1:L*L]],
        adsorptions=:full)

    cluster_triangular(nx, ny) = MLattice{1,TriangularLattice}(
        lattice_constant=1.0,
        supercell_dimensions=(nx, ny, 1),
        periodicity=(true, true, false),
        cutoff_radii=[1.1],
        components=:equal,
        adsorptions=:full)

    # ================================================================
    @testset "motif_distances" begin
        @test motif_distances([(0, 0), (1, 0), (0, 1)]) ≈ [1.0, 1.0, sqrt(2)]
        @test motif_distances([(0, 0), (1, 0), (0, 1), (1, 1)]) ≈
              [1.0, 1.0, 1.0, 1.0, sqrt(2), sqrt(2)]
        @test motif_distances([(0.0, 0.0, 0.0), (0.0, 0.0, 2.0)]) ≈ [2.0]
        @test_throws ArgumentError motif_distances([(0, 0)])
    end

    # ================================================================
    @testset "embedding counts on faithful cells" begin
        sq8 = cluster_square(8)

        # Right isoceles trio (1, 1, √2): four per plaquette
        t1 = enumerate_motif_embeddings(sq8, motif_distances([(0, 0), (1, 0), (0, 1)]);
                                        expected_count=256)
        @test length(t1) == 256
        counts = zeros(Int, 64)
        for e in t1, s in e
            counts[s] += 1
        end
        @test all(==(12), counts)   # per-site membership is uniform
        @test all(e -> e[1] < e[2] < e[3], t1)   # canonical form

        # Linear trio (1, 1, 2): two orientations per site pair
        @test length(enumerate_motif_embeddings(sq8, [1.0, 1.0, 2.0];
                                                expected_count=128)) == 128

        # Unit-square quattro via the template method: one per plaquette
        @test length(enumerate_motif_embeddings(sq8,
            [(0, 0), (1, 0), (0, 1), (1, 1)]; expected_count=64)) == 64
        # The multiset method warns for K ≥ 4 (homometric figures share
        # multisets); the unit square has no alias, so the counts agree
        q_ms = @test_logs (:warn, r"[Hh]omometric") match_mode = :any enumerate_motif_embeddings(
            sq8, motif_distances([(0, 0), (1, 0), (0, 1), (1, 1)]))
        @test length(q_ms) == 64

        # Pair signature (K = 2) reproduces the nearest-neighbor shell as
        # unordered pairs: 2M on the square lattice
        @test length(enumerate_motif_embeddings(sq8, [1.0];
                                                expected_count=128)) == 128

        # Triangular lattice, faithful cell: faces and linear trios
        tri44 = cluster_triangular(4, 4)   # M = 32
        @test length(enumerate_motif_embeddings(tri44, [1.0, 1.0, 1.0];
                                                expected_count=64)) == 64
        # (1, 1, 2) has K·d_max = 6 > C_x = 4: the wrap warning fires, and
        # the torus-convention count still equals 3M (ring-of-four subsets
        # coincide with the consecutive triples)
        tri_lin = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            tri44, [1.0, 1.0, 2.0])
        @test length(tri_lin) == 96

        # expected_count mismatches throw
        @test_throws ArgumentError enumerate_motif_embeddings(
            sq8, [1.0, 1.0, sqrt(2)]; expected_count=999)
        # Signature-length and value validation
        @test_throws ArgumentError enumerate_motif_embeddings(sq8, [1.0, 1.0])
        @test_throws ArgumentError enumerate_motif_embeddings(sq8, [1.0, 1.0, -2.0])
        @test_throws ArgumentError enumerate_motif_embeddings(sq8, [1.0, 1.0, 2.0]; tol=0.0)
    end

    # ================================================================
    @testset "wrap-around pathology (torus convention)" begin
        # The 18-site triangular cell: 36 faces plus 6 winding three-cycles,
        # the Study V5 diagnostic number, with the wrap warning
        tri33 = cluster_triangular(3, 3)   # M = 18
        t = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            tri33, [1.0, 1.0, 1.0])
        @test length(t) == 42

        # Square 4×4 with the linear trio: distance 2 = C/2 ties; the count
        # equals every 3-subset of each 4-ring, which coincides with 2M here
        sq4 = cluster_square(4)
        t4 = @test_logs (:warn, r"faithful quotient") match_mode = :any enumerate_motif_embeddings(
            sq4, [1.0, 1.0, 2.0])
        @test length(t4) == 32
    end

    # ================================================================
    @testset "homometric figures and tolerance semantics" begin
        # A homometric pair: non-congruent quattros sharing one distance
        # multiset. The template method separates the two families (each has
        # a mirror stabilizer of order 2, hence 4M embeddings); the multiset
        # method warns and enumerates both together.
        sq10 = cluster_square(10)
        A = [(0, 0), (1, 0), (0, 1), (2, 2)]
        B = [(0, 0), (1, 0), (2, 1), (2, 2)]
        @test motif_distances(A) ≈ motif_distances(B)
        eA = enumerate_motif_embeddings(sq10, A; expected_count=400)
        eB = enumerate_motif_embeddings(sq10, B; expected_count=400)
        @test isempty(intersect(Set(eA), Set(eB)))
        both = @test_logs (:warn, r"[Hh]omometric") match_mode = :any enumerate_motif_embeddings(
            sq10, motif_distances(A))
        @test length(both) == 800
        @test Set(both) == union(Set(eA), Set(eB))

        # Candidate pruning is a relaxation of acceptance: a distance within
        # tol of a signature entry that was merged into a nearby uniq
        # representative must still be found (regression: pruning at tol
        # from the representative silently dropped it)
        gen = MLattice{1,GenericLattice}(
            [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 1.0],
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.1, 0.0, 0.0)],
            (1, 1, 1),
            (false, false, false),
            [2.2],
            [[false, false, false]],
            [true, true, true])
        found = enumerate_motif_embeddings(gen, [1.0, 1.05, 2.1]; tol=0.06)
        @test length(found) == 1   # sides (1.0, 1.1, 2.1) accepted elementwise
    end

    # ================================================================
    @testset "cluster energy evaluation" begin
        Random.seed!(139)
        L = 6
        M = L * L
        sq6 = cluster_square(L)
        t1_embs = enumerate_motif_embeddings(sq6,
            motif_distances([(0, 0), (1, 0), (0, 1)]); expected_count=4M)
        q_embs = enumerate_motif_embeddings(sq6,
            [(0, 0), (1, 0), (0, 1), (1, 1)]; expected_count=M)
        V0, J1, J2 = -0.1, 0.05, 0.02
        Jt, Jq = 0.168, -0.120
        pair = GenericLatticeHamiltonian(V0, [J1, J2], u"eV")
        h = ClusterLatticeHamiltonian(pair,
            [ClusterInteraction(Jt * u"eV", t1_embs),
             ClusterInteraction(Jq * u"eV", q_embs)])

        # Empty lattice, single site, and full lattice (closed forms: every
        # site has 4 first- and 4 second-shell neighbors, 4M trio and M
        # quattro embeddings are all occupied)
        sq6.components[1] .= false
        @test interacting_energy(sq6, h) == 0.0u"eV"
        sq6.components[1][1] = true
        @test interacting_energy(sq6, h) ≈ V0 * u"eV" atol = 1e-14 * u"eV"
        sq6.components[1] .= true
        E_full = M * (V0 + 2J1 + 2J2) + 4M * Jt + M * Jq
        @test interacting_energy(sq6, h) ≈ E_full * u"eV" atol = 1e-10 * u"eV"

        # 50 random configurations against an independent brute force that
        # classifies pairs, trios, and quattros by integer squared
        # minimum-image distances; no code shared with the enumerator
        coord(s) = ((s - 1) % L, (s - 1) ÷ L)
        function d2(a, b)
            (ia, ja) = coord(a)
            (ib, jb) = coord(b)
            dx = min(mod(ia - ib, L), mod(ib - ia, L))
            dy = min(mod(ja - jb, L), mod(jb - ja, L))
            return dx^2 + dy^2
        end
        function brute_E(occ)
            s = findall(occ)
            n = length(s)
            E = V0 * n
            for x in 1:n, y in (x+1):n
                r2 = d2(s[x], s[y])
                r2 == 1 && (E += J1)
                r2 == 2 && (E += J2)
            end
            for x in 1:n, y in (x+1):n, z in (y+1):n
                r2s = sort([d2(s[x], s[y]), d2(s[x], s[z]), d2(s[y], s[z])])
                r2s == [1, 1, 2] && (E += Jt)
            end
            for w in 1:n, x in (w+1):n, y in (x+1):n, z in (y+1):n
                q = (s[w], s[x], s[y], s[z])
                r2s = sort([d2(q[a], q[b]) for a in 1:4 for b in (a+1):4])
                r2s == [1, 1, 1, 1, 2, 2] && (E += Jq)
            end
            return E
        end
        for _ in 1:50
            occ = rand(M) .< 0.5
            sq6.components[1] .= occ
            @test interacting_energy(sq6, h).val ≈ brute_E(occ) atol = 1e-10
        end

        # Translation and point-group images leave the energy unchanged
        occ0 = rand(M) .< 0.5
        translate(occ) = [occ[((i0 + 1) % L) + j0 * L + 1]
                          for j0 in 0:L-1 for i0 in 0:L-1]
        rotate(occ) = [occ[j0 + ((L - i0) % L) * L + 1]
                       for j0 in 0:L-1 for i0 in 0:L-1]
        sq6.components[1] .= occ0
        E0 = interacting_energy(sq6, h).val
        sq6.components[1] .= translate(occ0)
        @test interacting_energy(sq6, h).val ≈ E0 atol = 1e-12
        sq6.components[1] .= rotate(occ0)
        @test interacting_energy(sq6, h).val ≈ E0 atol = 1e-12

        # Fixed-N exact enumeration with a trio term: full N = 3 spectrum on
        # 4×4 against a hand-rolled triple loop
        lat4 = cluster_square(4)
        t1_4 = enumerate_motif_embeddings(lat4,
            motif_distances([(0, 0), (1, 0), (0, 1)]); expected_count=64)
        h4 = ClusterLatticeHamiltonian(
            GenericLatticeHamiltonian(0.0, [0.05, 0.0], u"eV"),
            [ClusterInteraction(0.168u"eV", t1_4)])
        lat4.components[1] .= false
        lat4.components[1][1:3] .= true
        df_exact, _ = exact_enumeration(lat4, h4)
        coord4(s) = ((s - 1) % 4, (s - 1) ÷ 4)
        function d24(a, b)
            (ia, ja) = coord4(a)
            (ib, jb) = coord4(b)
            dx = min(mod(ia - ib, 4), mod(ib - ia, 4))
            dy = min(mod(ja - jb, 4), mod(jb - ja, 4))
            return dx^2 + dy^2
        end
        brute3 = Dict{NTuple{3,Int},Float64}()
        for a in 1:16, b in (a+1):16, c in (b+1):16
            r2s = sort([d24(a, b), d24(a, c), d24(b, c)])
            brute3[(a, b, c)] = 0.05 * count(==(1), r2s) +
                                (r2s == [1, 1, 2] ? 0.168 : 0.0)
        end
        @test nrow(df_exact) == 560
        # Configuration-resolved: each enumerated configuration's energy
        # against the brute force keyed by its occupied sites (a multiset
        # comparison could not detect energies permuted among configurations)
        for row in eachrow(df_exact)
            occ_sites = NTuple{3,Int}(findall(row.config[1]))
            @test abs(ustrip(u"eV", row.energy) - brute3[occ_sites]) < 1e-10
        end

        # A Hamiltonian built for a larger lattice fails fast at liveset
        # construction instead of raising a BoundsError mid-run — and at the
        # raw-lattice sampler entry points, which never build a liveset
        w4 = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
        @test_throws ArgumentError LatticeGasWalkers(w4, h)
        wl_params = WangLandauParameters(energy_min=-8.0, energy_max=8.0,
                                         num_energy_bins=30, num_steps=10,
                                         max_iter=2, random_seed=7)
        @test_throws ArgumentError wang_landau(deepcopy(lat4), h, wl_params)
        @test_throws ArgumentError nvt_monte_carlo(MCNewSample(), deepcopy(lat4),
                                                   h, 300.0, 5, 7)

        # The converse mismatch — a Hamiltonian enumerated on a smaller
        # lattice, all of whose indices exist on the bigger one — is
        # undetectable from indices alone: construction succeeds and the
        # embeddings are geometric nonsense. Documented caller
        # responsibility: always enumerate on the lattice being sampled.
        w6 = [LatticeWalker(deepcopy(sq6), energy=0.0u"eV", iter=0)]
        @test LatticeGasWalkers(w6, h4) isa LatticeGasWalkers
    end

    # ================================================================
    @testset "Zhang O-Pd(100) m = 9 expansion: counts, adlayers, sampling" begin
        # Figure geometry transcribed from Zhang, Blum & Reuter, PRB 75,
        # 235406 (2007), Fig. 1 (arXiv:cond-mat/0701549). In the paper's own
        # labeling, V^t_1 is the LINEAR (1, 1, 2) trio and V^t_2 the
        # (1, 1, √2) right triangle; t3 is the scalene (1, √2, √5), t6 the
        # diagonal linear (√2, √2, 2√2), and q2 the "hut"
        # (0,0)-(1,0)-(2,0)+(1,1). FreeBird couplings are the negated GGA
        # values of their Table III (positive = repulsive).
        zhang_motifs = (
            t1=[(0, 0), (1, 0), (2, 0)],
            t2=[(0, 0), (1, 0), (0, 1)],
            t3=[(0, 1), (1, 0), (2, 0)],
            t6=[(0, 0), (1, 1), (2, 2)],
            q2=[(0, 0), (1, 0), (2, 0), (1, 1)])
        zhang_J = (t1=0.168, t2=-0.060, t3=0.048, t6=0.051, q2=-0.120)
        # Embeddings per site under standard full-orbit counting: 8 divided
        # by the order of the motif's point-group stabilizer
        zhang_per_site = (t1=2, t2=4, t3=8, t6=2, q2=4)

        function zhang_hamiltonian(L)
            lat = cluster_square(L; cutoffs=[1.1, 1.5, 2.1, 2.3])
            clusters = [ClusterInteraction(zhang_J[k] * u"eV",
                            enumerate_motif_embeddings(lat, zhang_motifs[k];
                                expected_count=zhang_per_site[k] * L * L))
                        for k in keys(zhang_motifs)]
            ham = ClusterLatticeHamiltonian(
                GenericLatticeHamiltonian(-1.249, [0.292, 0.090, -0.050, -0.010], u"eV"),
                clusters)
            return lat, ham
        end

        L12 = 12
        M12 = L12 * L12
        sq12, zhang = zhang_hamiltonian(L12)

        # Ordered adlayers: closed forms from hand-derived per-adatom motif
        # multiplicities (pins transcription, sign convention, and counting
        # simultaneously)
        set_phase!(pred) = (sq12.components[1] .=
            [pred((s - 1) % L12, (s - 1) ÷ L12) for s in 1:M12])

        # (3×3), 16 adatoms, no neighbor within √5: on-site only
        set_phase!((i, j) -> i % 3 == 0 && j % 3 == 0)
        @test interacting_energy(sq12, zhang).val ≈ 16 * (-1.249) atol = 1e-10
        # p(2×2), 36 adatoms: 2 third-shell pairs per adatom, no multi-body
        set_phase!((i, j) -> i % 2 == 0 && j % 2 == 0)
        @test interacting_energy(sq12, zhang).val ≈
              36 * (-1.249 + 2 * (-0.050)) atol = 1e-10
        # c(2×2), 72 adatoms: 2 second- and 2 third-shell pairs plus 2
        # diagonal-linear t6 trios per adatom
        set_phase!((i, j) -> iseven(i + j))
        @test interacting_energy(sq12, zhang).val ≈
              72 * (-1.249 + 2 * 0.090 + 2 * (-0.050) + 2 * 0.051) atol = 1e-10
        # p(2×1) stripes (every other column), 72 adatoms: pins J1 and the
        # linear t1 separately from t2 (per adatom: 1 vertical first-shell
        # pair, 2 third-shell pairs, 2 fourth-shell pairs, 1 vertical t1
        # trio; no √2 pairs, so t2, t3, t6, and q2 are all inactive)
        set_phase!((i, j) -> i % 2 == 0)
        @test interacting_energy(sq12, zhang).val ≈
              72 * (-1.249 + 0.292 + 2 * (-0.050) + 2 * (-0.010) + 0.168) atol = 1e-10
        # (1×1), 144 adatoms: every figure at its per-site multiplicity
        sq12.components[1] .= true
        @test interacting_energy(sq12, zhang).val ≈
              144 * (-1.249 +
                     2 * 0.292 + 2 * 0.090 + 2 * (-0.050) + 4 * (-0.010) +
                     2 * 0.168 + 4 * (-0.060) + 8 * 0.048 + 2 * 0.051 +
                     4 * (-0.120)) atol = 1e-10

        # The Θ ≤ 1/2 closed forms also match the paper's Table I DFT
        # binding energies to ≤ 2 meV per adatom under this counting
        # convention. The dense phases ((2×2)-3O, (1×1)) follow a different
        # multiplicity convention in the paper and are deliberately not
        # asserted against Table I.
        for (pred, n_ad, Eb_table) in (
                ((i, j) -> i % 3 == 0 && j % 3 == 0, 16, 1.249),
                ((i, j) -> i % 2 == 0 && j % 2 == 0, 36, 1.348),
                ((i, j) -> iseven(i + j), 72, 1.069))
            set_phase!(pred)
            @test isapprox(-interacting_energy(sq12, zhang).val / n_ad, Eb_table;
                           atol=0.003)
        end

        # The full Zhang Hamiltonian runs in GC-NS unmodified: seeded igref
        # run on 6×6, live-set energies match direct recomputation
        Random.seed!(1390)
        sq6z, z6 = zhang_hamiltonian(6)
        walkers = [LatticeWalker(deepcopy(sq6z), energy=0.0u"eV", iter=0)
                   for _ in 1:20]
        ls = LatticeGasWalkers(walkers, z6; assign_energy=false)
        params = IdealGasReferencedGCNSParameters(
            mc_steps=20, reference_fugacity=1.0, energy_perturbation=1e-9,
            allowed_fail_count=100_000)
        save = SaveEveryN("t_zh.csv", "t_zh.traj", "t_zh.ls",
                          1000000, 1000000, 1000000)
        df, fls, _ = ideal_gas_referenced_nested_sampling(
            ls, params, Int64(150),
            MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3), save)
        rm.(["t_zh.csv", "t_zh.traj", "t_zh.ls"], force=true)
        @test nrow(df) > 0
        for w in fls.walkers
            @test isapprox(w.energy.val,
                           interacting_energy(w.configuration, z6).val;
                           atol=1e-8)
        end
    end
end
