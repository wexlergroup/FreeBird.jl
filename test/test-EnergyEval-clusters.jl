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
end
