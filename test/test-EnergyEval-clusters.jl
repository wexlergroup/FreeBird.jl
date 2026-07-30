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

        # Unit-square quattro: one per plaquette
        @test length(enumerate_motif_embeddings(sq8,
            motif_distances([(0, 0), (1, 0), (0, 1), (1, 1)]);
            expected_count=64)) == 64

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
    @testset "cluster energy evaluation" begin
        Random.seed!(139)
        L = 6
        M = L * L
        sq6 = cluster_square(L)
        t1_embs = enumerate_motif_embeddings(sq6,
            motif_distances([(0, 0), (1, 0), (0, 1)]); expected_count=4M)
        q_embs = enumerate_motif_embeddings(sq6,
            motif_distances([(0, 0), (1, 0), (0, 1), (1, 1)]); expected_count=M)
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
        brute3 = Float64[]
        for a in 1:16, b in (a+1):16, c in (b+1):16
            r2s = sort([d24(a, b), d24(a, c), d24(b, c)])
            push!(brute3, 0.05 * count(==(1), r2s) +
                          (r2s == [1, 1, 2] ? 0.168 : 0.0))
        end
        lib3 = sort([ustrip(u"eV", e) for e in df_exact.energy])
        @test length(lib3) == 560
        @test maximum(abs.(lib3 .- sort(brute3))) < 1e-10

        # A Hamiltonian built for a larger lattice fails fast at liveset
        # construction instead of raising a BoundsError mid-run
        w4 = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
        @test_throws ArgumentError LatticeGasWalkers(w4, h)
    end
end
