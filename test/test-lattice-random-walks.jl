@testset "Monte Carlo Moves tests" begin
    @testset "lattice random walk tests" begin

        sl = SLattice{SquareLattice}(components=[[1,2,3,4]])
        new_sl = lattice_random_walk!(deepcopy(sl))

        @test occupied_site_count(sl) == occupied_site_count(new_sl)
        @test length(sl.components[1]) == length(new_sl.components[1])

        ml = MLattice{2,SquareLattice}()
        new_ml = lattice_random_walk!(deepcopy(ml))

        @test occupied_site_count(ml) == occupied_site_count(new_ml)
        @test length(ml.components[1]) == length(new_ml.components[1])
        @test length(ml.components[2]) == length(new_ml.components[2])
    end

    @testset "Swap occupied sites across components tests" begin
        ml = MLattice{2,SquareLattice}(components=[[1,3],[2,4]])
        new_ml = deepcopy(ml)
        MonteCarloMoves.swap_occupied_sites_across_components!(new_ml, 1, 2)

        @test new_ml.components[1][1] == 0
        @test new_ml.components[1][2] == 1
        @test new_ml.components[2][1] == 1
        @test new_ml.components[2][2] == 0
    end

    @testset "Swap empty and occupied sites tests" begin
        ml = MLattice{2,SquareLattice}(components=[[1],[2]])
        new_ml = deepcopy(ml)
        MonteCarloMoves.swap_empty_occupied_sites!(new_ml, 1, 2)

        @test new_ml.components[1][1] == 0
        @test new_ml.components[1][2] == 1
        @test new_ml.components[2][1] == 1
        @test new_ml.components[2][2] == 0
    end

    @testset "geometric cluster move tests" begin

        @testset "reflection map is self-inverse (2D)" begin
            Lx, Ly, Lz = 6, 6, 1
            for pivot_gx in [0, 1, 3, 5], pivot_gy in [0, 2, 4, 5]
                for site in 1:(Lx * Ly * Lz)
                    r = MonteCarloMoves._reflect_site(site, pivot_gx, pivot_gy, Lx, Ly, Lz)
                    rr = MonteCarloMoves._reflect_site(r, pivot_gx, pivot_gy, Lx, Ly, Lz)
                    @test rr == site
                end
            end
        end

        @testset "reflection map is self-inverse (3D)" begin
            Lx, Ly, Lz = 4, 4, 3
            for pivot_gx in [0, 1, 3], pivot_gy in [0, 2, 3]
                for site in 1:(Lx * Ly * Lz)
                    r = MonteCarloMoves._reflect_site(site, pivot_gx, pivot_gy, Lx, Ly, Lz)
                    rr = MonteCarloMoves._reflect_site(r, pivot_gx, pivot_gy, Lx, Ly, Lz)
                    @test rr == site
                end
            end
        end

        @testset "reflection preserves z-coordinate (3D)" begin
            Lx, Ly, Lz = 4, 4, 3
            for pivot_gx in [0, 1, 2], pivot_gy in [0, 1, 3]
                for site in 1:(Lx * Ly * Lz)
                    _, _, gz_orig = MonteCarloMoves._site_to_grid(site, Lx, Ly)
                    r = MonteCarloMoves._reflect_site(site, pivot_gx, pivot_gy, Lx, Ly, Lz)
                    _, _, gz_refl = MonteCarloMoves._site_to_grid(r, Lx, Ly)
                    @test gz_refl == gz_orig
                end
            end
        end

        @testset "reflection wraps correctly under PBC" begin
            Lx, Ly, Lz = 4, 4, 1
            # Pivot at (0,0): R(1,0) -> (2*0 - 1 mod 4, 0) = (3, 0)
            site_10 = MonteCarloMoves._grid_to_site(1, 0, 0, Lx, Ly)  # site at grid (1,0,0)
            reflected = MonteCarloMoves._reflect_site(site_10, 0, 0, Lx, Ly, Lz)
            gx, gy, gz = MonteCarloMoves._site_to_grid(reflected, Lx, Ly)
            @test gx == 3
            @test gy == 0
            @test gz == 0

            # Pivot at (2,2): R(0,0) -> (4 mod 4, 4 mod 4) = (0, 0) — fixed point
            site_00 = MonteCarloMoves._grid_to_site(0, 0, 0, Lx, Ly)
            reflected = MonteCarloMoves._reflect_site(site_00, 2, 2, Lx, Ly, Lz)
            gx, gy, gz = MonteCarloMoves._site_to_grid(reflected, Lx, Ly)
            @test gx == 0
            @test gy == 0

            # Pivot at (1,1): R(3,3) -> (2-3 mod 4, 2-3 mod 4) = (3, 3) — check wrap
            site_33 = MonteCarloMoves._grid_to_site(3, 3, 0, Lx, Ly)
            reflected = MonteCarloMoves._reflect_site(site_33, 1, 1, Lx, Ly, Lz)
            gx, gy, gz = MonteCarloMoves._site_to_grid(reflected, Lx, Ly)
            @test gx == mod(2*1 - 3, 4)  # 3
            @test gy == mod(2*1 - 3, 4)  # 3
        end

        @testset "site ↔ grid round-trip (2D)" begin
            Lx, Ly, Lz = 5, 7, 1
            for site in 1:(Lx * Ly * Lz)
                gx, gy, gz = MonteCarloMoves._site_to_grid(site, Lx, Ly)
                @test MonteCarloMoves._grid_to_site(gx, gy, gz, Lx, Ly) == site
                @test 0 <= gx < Lx
                @test 0 <= gy < Ly
                @test gz == 0
            end
        end

        @testset "site ↔ grid round-trip (3D)" begin
            Lx, Ly, Lz = 4, 4, 3
            for site in 1:(Lx * Ly * Lz)
                gx, gy, gz = MonteCarloMoves._site_to_grid(site, Lx, Ly)
                @test MonteCarloMoves._grid_to_site(gx, gy, gz, Lx, Ly) == site
                @test 0 <= gx < Lx
                @test 0 <= gy < Ly
                @test 0 <= gz < Lz
            end
        end

        @testset "particle count preserved (C=1)" begin
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(8, 8, 1),
                components=[[1, 2, 3, 10, 15, 20, 30, 40, 50, 60]]
            )
            original_count = occupied_site_count(sl)
            for _ in 1:50
                geometric_cluster_swap!(sl, 0.3)
                @test occupied_site_count(sl) == original_count
            end
        end

        @testset "particle count preserved (C=2)" begin
            ml = MLattice{2,SquareLattice}(
                supercell_dimensions=(6, 6, 1),
                components=[[1, 3, 5, 7, 9, 11], [2, 4, 6, 8, 10, 12]]
            )
            original_counts = occupied_site_count(ml)
            for _ in 1:50
                geometric_cluster_swap!(ml, 0.4)
                @test occupied_site_count(ml) == original_counts
            end
        end

        @testset "self-inverse for fixed cluster" begin
            using Random
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(8, 8, 1),
                components=[[1, 5, 10, 20, 30, 40, 50, 60]]
            )
            original_components = deepcopy(sl.components)

            # Apply with seeded RNG, then apply again with same seed
            seed = 42
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)

            @test sl.components == original_components
        end

        @testset "self-inverse for fixed cluster (C=2)" begin
            using Random
            ml = MLattice{2,SquareLattice}(
                supercell_dimensions=(6, 6, 1),
                components=[[1, 3, 5, 7, 9], [2, 4, 6, 8, 10]]
            )
            original_components = deepcopy(ml.components)

            seed = 123
            Random.seed!(seed)
            geometric_cluster_swap!(ml, 0.4)
            Random.seed!(seed)
            geometric_cluster_swap!(ml, 0.4)

            @test ml.components == original_components
        end

        @testset "particle count preserved 3D (C=1)" begin
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(4, 4, 3),
                periodicity=(true, true, false),
                components=[[1, 5, 20, 35]]
            )
            original_count = occupied_site_count(sl)
            for _ in 1:50
                geometric_cluster_swap!(sl, 0.3)
                @test occupied_site_count(sl) == original_count
            end
        end

        @testset "self-inverse for fixed cluster 3D" begin
            using Random
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(4, 4, 3),
                periodicity=(true, true, false),
                components=[[1, 5, 17, 33, 40]]
            )
            original_components = deepcopy(sl.components)

            seed = 42
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)

            @test sl.components == original_components
        end

        @testset "cluster move can change configuration" begin
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(8, 8, 1),
                components=[[1, 2, 3, 4, 5, 6, 7, 8]]
            )
            original_components = deepcopy(sl.components)
            changed = false
            for _ in 1:100
                test_sl = deepcopy(sl)
                geometric_cluster_swap!(test_sl, 0.5)
                if test_sl.components != original_components
                    changed = true
                    break
                end
            end
            @test changed
        end

        @testset "multi-site basis guard" begin
            using Random
            # The keyword constructor admits an arbitrary basis, so the
            # index-space reflection's single-basis contract gets the same
            # loud guard the triangular method has. (Basis offset and cutoff
            # chosen so the intra-basis shell is unambiguous — a warning-free
            # construction.)
            two_site = MLattice{1,SquareLattice}(
                basis=[(0.0, 0.0, 0.0), (0.25, 0.25, 0.0)],
                supercell_dimensions=(4, 4, 1),
                cutoff_radii=[0.5],
                components=[[1, 2, 3]]
            )
            @test_throws ArgumentError geometric_cluster_swap!(two_site, 0.3)

            # The guard consumes no randomness: a guarded throw leaves the
            # global RNG stream untouched.
            Random.seed!(7)
            probe = rand(UInt64)
            Random.seed!(7)
            try
                geometric_cluster_swap!(two_site, 0.3)
            catch
            end
            @test rand(UInt64) === probe

            # On a single-site-basis lattice the guarded method's same-seed
            # trajectory is identical to the pre-guard body, replicated here
            # through the same internal helpers.
            sl = SLattice{SquareLattice}(
                supercell_dimensions=(8, 8, 1),
                components=[[1, 2, 3, 10, 15, 20, 30, 40, 50, 60]]
            )
            replica = deepcopy(sl)
            Random.seed!(4242)
            for _ in 1:25
                geometric_cluster_swap!(sl, 0.3)
            end
            Random.seed!(4242)
            for _ in 1:25
                Lx, Ly, Lz = replica.supercell_dimensions
                pivot_gx = rand(0:Lx-1)
                pivot_gy = rand(0:Ly-1)
                seed_site = rand(1:num_sites(replica))
                reflect = s -> MonteCarloMoves._reflect_site(s, pivot_gx, pivot_gy, Lx, Ly, Lz)
                cluster = MonteCarloMoves._build_geometric_cluster(replica, seed_site, reflect, 0.3)
                for (a, b) in cluster
                    if a != b
                        replica.components[1][a], replica.components[1][b] =
                            replica.components[1][b], replica.components[1][a]
                    end
                end
            end
            @test sl.components == replica.components
        end

        @testset "GenericLattice guard" begin
            using Random
            # Geometric cluster moves carry square and triangular reflection
            # maps only; a GenericLattice configuration gets the same loud
            # guard so a cluster-armed routine fails at the move's definition
            # with a descriptive ArgumentError instead of a raw MethodError
            # from inside the sampling loop (the keyword MCMixedMoves
            # constructor arms cluster moves by default).
            gen = MLattice{1,GenericLattice}(
                [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
                [(0.0, 0.0, 0.0)],
                (3, 1, 1),
                (false, false, false),
                [1.1],
                [[true, false, false]],
                fill(false, 3))
            err = try
                geometric_cluster_swap!(gen, 0.3)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("square and triangular", err.msg)
            @test occursin("GenericLattice", err.msg)

            # The guard consumes no randomness: a guarded throw leaves the
            # global RNG stream untouched.
            Random.seed!(11)
            probe = rand(UInt64)
            Random.seed!(11)
            try
                geometric_cluster_swap!(gen, 0.3)
            catch
            end
            @test rand(UInt64) === probe

            # End to end: the cluster-armed keyword routine surfaces the
            # guard's ArgumentError through the sampling step, while the
            # positional back-compat form stays clusters-free and completes.
            ham = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
            walkers = [LatticeWalker(deepcopy(gen), energy=0.0u"eV", iter=0) for _ in 1:4]
            ls = LatticeGasWalkers(walkers, ham)
            params = NestedSamplingParameters(mc_steps=10)
            @test_throws ArgumentError SamplingSchemes.nested_sampling_step!(ls, params, MCMixedMoves())
            Random.seed!(99)
            walkers2 = [LatticeWalker(deepcopy(gen), energy=0.0u"eV", iter=0) for _ in 1:4]
            ls2 = LatticeGasWalkers(walkers2, ham)
            save_strategy = SaveEveryN(n_traj=10^6, n_snap=10^6, n_info=10^6)
            df, ls3, _ = nested_sampling(ls2, NestedSamplingParameters(mc_steps=10), 30, MCMixedMoves(5, 1), save_strategy)
            @test size(df, 1) >= 1
            @test length(ls3.walkers) == 4
        end

        # ---- triangular lattice (two-site centered-rectangular basis) ----

        @testset "triangular reflection: half-grid round-trip" begin
            for (nx, ny) in ((3, 3), (6, 4))
                M = 2 * nx * ny
                for s in 1:M
                    hx, hy = MonteCarloMoves._tri_site_to_halfgrid(s, nx)
                    # The centered-rectangular condition: hx ≡ hy (mod 2)
                    @test (hx & 1) == (hy & 1)
                    @test 0 <= hx < 2 * nx
                    @test 0 <= hy < 2 * ny
                    @test MonteCarloMoves._tri_halfgrid_to_site(hx, hy, nx) == s
                end
            end
        end

        @testset "triangular reflection: involution and bijection" begin
            using Random
            Random.seed!(5)
            for (nx, ny) in ((3, 3), (6, 4))
                M = 2 * nx * ny
                for _ in 1:25
                    h1 = MonteCarloMoves._tri_site_to_halfgrid(rand(1:M), nx)
                    h2 = MonteCarloMoves._tri_site_to_halfgrid(rand(1:M), nx)
                    hpx, hpy = h1[1] + h2[1], h1[2] + h2[2]
                    σ = [MonteCarloMoves._tri_reflect_site(s, hpx, hpy, nx, ny)
                         for s in 1:M]
                    @test sort(σ) == collect(1:M)          # bijection on sites
                    @test all(σ[σ[s]] == s for s in 1:M)   # involution
                end
            end
        end

        @testset "triangular reflection is a lattice symmetry" begin
            # Geometric check against the stored positions, independent of
            # the half-grid index arithmetic: 2·pivot − r(s) − r(σ(s)) must
            # be a supercell lattice vector (periods Ax = nx·a, Ay = ny·√3·a)
            using Random
            Random.seed!(6)
            for (nx, ny) in ((3, 3), (6, 4))
                lat = MLattice{1,TriangularLattice}(
                    supercell_dimensions=(nx, ny, 1),
                    cutoff_radii=[1.1],
                    components=[[1]],
                    adsorptions=:full)
                M = 2 * nx * ny
                Ax = nx * 1.0
                Ay = ny * sqrt(3)
                for _ in 1:10
                    s1 = rand(1:M)
                    s2 = rand(1:M)
                    h1 = MonteCarloMoves._tri_site_to_halfgrid(s1, nx)
                    h2 = MonteCarloMoves._tri_site_to_halfgrid(s2, nx)
                    px = (lat.positions[s1, 1] + lat.positions[s2, 1]) / 2
                    py = (lat.positions[s1, 2] + lat.positions[s2, 2]) / 2
                    for s in 1:M
                        r = MonteCarloMoves._tri_reflect_site(
                            s, h1[1] + h2[1], h1[2] + h2[2], nx, ny)
                        fx = (2 * px - lat.positions[s, 1] - lat.positions[r, 1]) / Ax
                        fy = (2 * py - lat.positions[s, 2] - lat.positions[r, 2]) / Ay
                        @test isapprox(fx, round(fx), atol=1e-9)
                        @test isapprox(fy, round(fy), atol=1e-9)
                    end
                end
            end
        end

        @testset "triangular basis-sublattice rule" begin
            nx, ny = 3, 3
            M = 2 * nx * ny
            basis_of(s) = (s - 1) % 2
            # Site pivots (s1 == s2, parity sum even): basis index preserved
            for s1 in (1, 4, 18)
                h = MonteCarloMoves._tri_site_to_halfgrid(s1, nx)
                for s in 1:M
                    r = MonteCarloMoves._tri_reflect_site(s, 2 * h[1], 2 * h[2], nx, ny)
                    @test basis_of(r) == basis_of(s)
                end
            end
            # Midpoint pivot with odd parity sum: basis 0 and basis 1 exchange
            h1 = MonteCarloMoves._tri_site_to_halfgrid(1, nx)   # basis 0
            h2 = MonteCarloMoves._tri_site_to_halfgrid(2, nx)   # basis 1
            @test isodd(h1[1] + h2[1])
            for s in 1:M
                r = MonteCarloMoves._tri_reflect_site(
                    s, h1[1] + h2[1], h1[2] + h2[2], nx, ny)
                @test basis_of(r) == 1 - basis_of(s)
            end
        end

        @testset "particle count preserved (triangular, C=1)" begin
            sl = SLattice{TriangularLattice}(
                supercell_dimensions=(6, 4, 1),
                cutoff_radii=[1.1],
                components=[[1, 2, 5, 10, 17, 24, 33, 40]],
                adsorptions=:full)
            original_count = occupied_site_count(sl)
            for _ in 1:50
                geometric_cluster_swap!(sl, 0.3)
                @test occupied_site_count(sl) == original_count
            end
        end

        @testset "particle count preserved (triangular, C=2)" begin
            ml = MLattice{2,TriangularLattice}(
                supercell_dimensions=(3, 3, 1),
                cutoff_radii=[1.1],
                components=[[1, 3, 5, 7], [2, 4, 6, 8]],
                adsorptions=:full)
            original_counts = occupied_site_count(ml)
            for _ in 1:50
                geometric_cluster_swap!(ml, 0.4)
                @test occupied_site_count(ml) == original_counts
            end
        end

        @testset "self-inverse for fixed cluster (triangular)" begin
            using Random
            sl = SLattice{TriangularLattice}(
                supercell_dimensions=(6, 4, 1),
                cutoff_radii=[1.1],
                components=[[1, 5, 10, 20, 30, 40]],
                adsorptions=:full)
            original_components = deepcopy(sl.components)

            seed = 42
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)
            Random.seed!(seed)
            geometric_cluster_swap!(sl, 0.3)

            @test sl.components == original_components
        end

        @testset "cluster move can change configuration (triangular)" begin
            sl = SLattice{TriangularLattice}(
                supercell_dimensions=(6, 4, 1),
                cutoff_radii=[1.1],
                components=[[1, 2, 3, 4, 5, 6]],
                adsorptions=:full)
            original_components = deepcopy(sl.components)
            changed = false
            for _ in 1:100
                test_sl = deepcopy(sl)
                geometric_cluster_swap!(test_sl, 0.5)
                if test_sl.components != original_components
                    changed = true
                    break
                end
            end
            @test changed
        end

        @testset "triangular guards" begin
            # One-site basis: the half-grid index arithmetic does not apply
            lat1b = MLattice{1,TriangularLattice}(
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(3, 3, 1),
                cutoff_radii=[1.1],
                components=[[1]],
                adsorptions=:full)
            @test_throws ArgumentError geometric_cluster_swap!(lat1b, 0.3)

            # Three-dimensional supercell: not supported on triangular
            lat3d = MLattice{1,TriangularLattice}(
                supercell_dimensions=(3, 3, 2),
                periodicity=(true, true, true),
                cutoff_radii=[1.1],
                components=[[1]],
                adsorptions=:full)
            @test_throws ArgumentError geometric_cluster_swap!(lat3d, 0.3)
        end

    end
end
