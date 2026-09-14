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

            # Stacked cell with a non-periodic third axis: the layer
            # reflection wraps on the c period, so this keeps throwing
            lat3d_open = MLattice{1,TriangularLattice}(
                supercell_dimensions=(3, 3, 2),
                periodicity=(true, true, false),
                cutoff_radii=[1.1],
                components=[[1]],
                adsorptions=:full)
            @test_throws ArgumentError geometric_cluster_swap!(lat3d_open, 0.3)

            # Disclosure: the assertion shipped here expected an ArgumentError
            # from this c-periodic (3, 3, 2) cell, stacked triangular cells
            # being unsupported. The three-dimensional point inversion accepts
            # it, so the assertion is replaced by the new contract: the move
            # runs and conserves the particle count.
            lat3d = MLattice{1,TriangularLattice}(
                supercell_dimensions=(3, 3, 2),
                periodicity=(true, true, true),
                cutoff_radii=[1.1],
                components=[[1, 7, 20]],
                adsorptions=:full)
            geometric_cluster_swap!(lat3d, 0.3)
            @test occupied_site_count(lat3d) == [3]
        end

        @testset "stacked triangular cluster moves" begin
            using Random
            using Unitful
            using DataFrames

            # Fixtures. Aligned stacks come from the keyword constructor with
            # `interlayer_spacing`; at the default cutoffs shell 2 holds the
            # axial interlayer neighbours (two at nz >= 3). The (6, 4, 2)
            # aligned cells are built with the in-plane shell only, since at
            # nz = 2 the two axial neighbours are the same site and the
            # default neighbour path would warn; those cells feed only the
            # map tests. Offset stacks, each layer shifted by the triangle
            # centre (the ABC stacking of close-packed layers), are built
            # through the inner constructor with the third lattice vector
            # (1/2, sqrt(3)/6, h) and `image_multiplicity=true` (the supercell
            # is non-orthogonal; see the geometric_cluster_swap! docstring);
            # at h = 1 and cutoffs [1.05, 1.2], shell 1 is the six in-plane
            # neighbours and shell 2 the six adjacent-layer neighbours
            # (verified by execution). Both stackings are c-periodic.
            B_of(nx, ny) = 2 * nx * ny
            layer_of(s, nx, ny) = (s - 1) ÷ B_of(nx, ny)
            basis_of(s) = (s - 1) % 2
            function aligned_tri(nx, ny, nz; h=1.2, cutoffs=[1.1, 1.25],
                                 comps=[[1]])
                MLattice{length(comps),TriangularLattice}(
                    supercell_dimensions=(nx, ny, nz),
                    periodicity=(true, true, true),
                    interlayer_spacing=h,
                    cutoff_radii=cutoffs,
                    components=comps,
                    adsorptions=:full)
            end
            function offset_tri(nx, ny, nz; h=1.0, cutoffs=[1.05, 1.2],
                                comps=[[1]])
                M = 2 * nx * ny * nz
                occ = [zeros(Bool, M) for _ in comps]
                for (c, sites) in enumerate(comps)
                    occ[c][sites] .= true
                end
                MLattice{length(comps),TriangularLattice}(
                    [1.0 0.0 0.5; 0.0 sqrt(3) sqrt(3)/6; 0.0 0.0 h],
                    [(0.0, 0.0, 0.0), (0.5, sqrt(3)/2, 0.0)],
                    (nx, ny, nz), (true, true, true), cutoffs,
                    occ, ones(Bool, M); image_multiplicity=true)
            end
            # The map under test, through the pivot midpoint of sites s1, s2
            function reflect3d(s, s1, s2, nx, ny, nz)
                h1 = MonteCarloMoves._tri_site_to_halfgrid(s1, nx, ny)
                h2 = MonteCarloMoves._tri_site_to_halfgrid(s2, nx, ny)
                return MonteCarloMoves._tri_reflect_site_3d(
                    s, h1[1] + h2[1], h1[2] + h2[2], h1[3] + h2[3], nx, ny, nz)
            end

            @testset "stacked reflection: involution and bijection" begin
                Random.seed!(7201)
                for lat in (aligned_tri(4, 4, 3), offset_tri(4, 4, 3),
                            aligned_tri(6, 4, 2; cutoffs=[1.1]), offset_tri(6, 4, 2))
                    nx, ny, nz = lat.supercell_dimensions
                    M = num_sites(lat)
                    for _ in 1:25
                        s1 = rand(1:M)
                        s2 = rand(1:M)
                        σ = [reflect3d(s, s1, s2, nx, ny, nz) for s in 1:M]
                        @test sort(σ) == collect(1:M)          # bijection on sites
                        @test all(σ[σ[s]] == s for s in 1:M)   # involution
                    end
                end
            end

            @testset "stacked reflection is a lattice symmetry" begin
                # Geometric check against the stored positions, independent of
                # the index arithmetic: 2·pivot − r(s) − r(σ(s)) must be a
                # supercell lattice vector (fractional coordinates integer),
                # on both stackings
                Random.seed!(7202)
                for lat in (aligned_tri(4, 4, 3), offset_tri(4, 4, 3),
                            aligned_tri(6, 4, 2; cutoffs=[1.1]), offset_tri(6, 4, 2))
                    nx, ny, nz = lat.supercell_dimensions
                    M = num_sites(lat)
                    A = lat.lattice_vectors * [nx 0 0; 0 ny 0; 0 0 nz]
                    Ainv = inv(A)
                    for _ in 1:10
                        s1 = rand(1:M)
                        s2 = rand(1:M)
                        pv = lat.positions[s1, :] .+ lat.positions[s2, :]
                        for s in 1:M
                            r = reflect3d(s, s1, s2, nx, ny, nz)
                            f = Ainv * (pv .- lat.positions[s, :] .- lat.positions[r, :])
                            @test all(isapprox.(f, round.(f), atol=1e-9))
                        end
                    end
                end
            end

            @testset "stacked reflection layer rule" begin
                # Any pivot sends layer k to kp − k (mod nz), an involution
                # on layers. A site pivot (kp = 2 k1) preserves the basis
                # index and fixes its own layer, sending layer k1 + d to
                # k1 − d; a midpoint pivot with an odd layer sum pairs the
                # layers whose indices sum to kp (mod nz), which at nz = 3
                # fixes the third layer and at nz = 2 is the strict exchange
                # 0 <-> 1
                for lat in (aligned_tri(4, 4, 3), offset_tri(4, 4, 3))
                    nx, ny, nz = lat.supercell_dimensions
                    M = num_sites(lat)
                    B = B_of(nx, ny)
                    for s1 in (1, B + 4, 2 * B + 18)
                        k1 = layer_of(s1, nx, ny)
                        for s in 1:M
                            r = reflect3d(s, s1, s1, nx, ny, nz)
                            @test basis_of(r) == basis_of(s)
                            @test layer_of(r, nx, ny) == mod(2 * k1 - layer_of(s, nx, ny), nz)
                        end
                        @test all(layer_of(reflect3d(s, s1, s1, nx, ny, nz), nx, ny) == k1
                                  for s in (k1 * B + 1):((k1 + 1) * B))
                    end
                    s1, s2 = 1, B + 2
                    kp = layer_of(s1, nx, ny) + layer_of(s2, nx, ny)
                    @test isodd(kp)
                    for s in 1:M
                        r = reflect3d(s, s1, s2, nx, ny, nz)
                        @test layer_of(r, nx, ny) == mod(kp - layer_of(s, nx, ny), nz)
                    end
                    # At nz = 3 the layer not paired by kp = 1 is fixed
                    @test all(layer_of(reflect3d(s, s1, s2, nx, ny, nz), nx, ny) == 2
                              for s in (2 * B + 1):(3 * B))
                end
                for lat in (aligned_tri(6, 4, 2; cutoffs=[1.1]), offset_tri(6, 4, 2))
                    nx, ny, nz = lat.supercell_dimensions
                    M = num_sites(lat)
                    B = B_of(nx, ny)
                    s1, s2 = 3, B + 5
                    @test isodd(layer_of(s1, nx, ny) + layer_of(s2, nx, ny))
                    for s in 1:M
                        r = reflect3d(s, s1, s2, nx, ny, nz)
                        @test layer_of(r, nx, ny) == 1 - layer_of(s, nx, ny)
                    end
                end
            end

            @testset "reduction to the single-layer map at nz = 1" begin
                nx, ny = 6, 4
                M = 2 * nx * ny
                @test all(MonteCarloMoves._tri_site_to_halfgrid(s, nx, ny) ==
                          (MonteCarloMoves._tri_site_to_halfgrid(s, nx)..., 0)
                          for s in 1:M)
                for s1 in 1:M, s2 in 1:M
                    h1 = MonteCarloMoves._tri_site_to_halfgrid(s1, nx)
                    h2 = MonteCarloMoves._tri_site_to_halfgrid(s2, nx)
                    hpx, hpy = h1[1] + h2[1], h1[2] + h2[2]
                    @test all(MonteCarloMoves._tri_reflect_site_3d(
                                  s, hpx, hpy, 0, nx, ny, 1) ==
                              MonteCarloMoves._tri_reflect_site(s, hpx, hpy, nx, ny)
                              for s in 1:M)
                end

                # Seeded single-layer cluster walk, captured on the
                # single-layer method shipped before this change (dev
                # 60113fd1) and reproduced identically across two Julia
                # processes. The ceiling binds: 36 of 60 proposals are
                # rejected at the first seed and 15 at the second, so the
                # revert path is inside the pinned stream. Every recorded
                # float comes from the fixed nested scalar accumulation of
                # the lattice energy plus one perturbation product per step,
                # with no vectorized reduction, so the exact pins are expected
                # to hold on every CI leg. Disclosed fallback: should a leg
                # falsify that on the energy digits, those two pins drop to
                # rtol 1e-12 while the integer and occupation pins stay exact.
                pin_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
                function pin_walk(seed)
                    sl = SLattice{TriangularLattice}(
                        supercell_dimensions=(6, 4, 1),
                        cutoff_radii=[1.1, 1.8],
                        components=[[1, 5, 10, 20, 30, 40, 13, 27]],
                        adsorptions=:full)
                    w = LatticeWalker(sl, energy=interacting_energy(sl, pin_ham), iter=0)
                    Random.seed!(seed)
                    acc, rate, w = MC_cluster_walk!(60, w, pin_ham, -0.36, 0.3;
                                                    energy_perturb=1e-9)
                    return acc, rate, w.energy.val,
                           findall(w.configuration.components[1])
                end
                acc1, rate1, e1, occ1 = pin_walk(7101)
                @test acc1 == true
                @test rate1 == 0.4
                @test e1 == -0.4225000004312682
                @test occ1 == [8, 9, 10, 11, 12, 14, 21, 42]
                acc2, rate2, e2, occ2 = pin_walk(7102)
                @test acc2 == true
                @test rate2 == 0.75
                @test e2 == -0.39000000014891933
                @test occ2 == [6, 9, 16, 19, 26, 27, 29, 30]
                # Same-process replay identity
                @test pin_walk(7101) == (acc1, rate1, e1, occ1)
            end

            @testset "particle count preserved (stacked triangular)" begin
                Random.seed!(7204)
                ham1 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.006], u"eV")
                sites1 = [1, 2, 5, 10, 17, 24, 33, 40, 50, 61, 77, 90]
                for lat in (aligned_tri(4, 4, 3; comps=[sites1]),
                            offset_tri(4, 4, 3; comps=[sites1]))
                    w = LatticeWalker(lat, energy=interacting_energy(lat, ham1), iter=0)
                    n0 = occupied_site_count(lat)
                    for _ in 1:25
                        MC_cluster_walk!(4, w, ham1, Inf, 0.3)
                        @test occupied_site_count(w.configuration) == n0
                    end
                end
                # Two components: every off-diagonal and diagonal coupling set,
                # the matrix built directly (no flattened-vector constructor)
                h11 = GenericLatticeHamiltonian(-0.04, [-0.01, -0.006], u"eV")
                h12 = GenericLatticeHamiltonian(-0.03, [-0.005, -0.002], u"eV")
                h22 = GenericLatticeHamiltonian(-0.02, [-0.008, -0.004], u"eV")
                ham2 = MLatticeHamiltonian{2,2,typeof(0.0u"eV")}(
                    reshape([h11, h12, h12, h22], 2, 2))
                sitesA = [1, 3, 5, 7, 34, 36, 66, 70]
                sitesB = [2, 4, 6, 8, 35, 37, 67, 71]
                for lat in (aligned_tri(4, 4, 3; comps=[sitesA, sitesB]),
                            offset_tri(4, 4, 3; comps=[sitesA, sitesB]))
                    w = LatticeWalker(lat, energy=interacting_energy(lat, ham2), iter=0)
                    n0 = occupied_site_count(lat)
                    for _ in 1:25
                        MC_cluster_walk!(4, w, ham2, Inf, 0.4)
                        @test occupied_site_count(w.configuration) == n0
                    end
                end
            end

            @testset "self-inverse and layer crossing (stacked triangular)" begin
                sites = [1, 5, 10, 20, 30, 40, 55, 70, 85]
                for lat in (aligned_tri(4, 4, 3; comps=[sites]),
                            offset_tri(4, 4, 3; comps=[sites]))
                    nx, ny, nz = lat.supercell_dimensions
                    original = deepcopy(lat.components)
                    Random.seed!(7205)
                    geometric_cluster_swap!(lat, 0.3)
                    Random.seed!(7205)
                    geometric_cluster_swap!(lat, 0.3)
                    @test lat.components == original
                    # The layer reflection is exercised: some applied pair
                    # joins sites in different layers
                    crossed = false
                    rec = Tuple{Int,Int}[]
                    for _ in 1:200
                        empty!(rec)
                        geometric_cluster_swap!(lat, 0.3; record=rec)
                        if any(layer_of(a, nx, ny) != layer_of(b, nx, ny) for (a, b) in rec)
                            crossed = true
                            break
                        end
                    end
                    @test crossed
                end
            end

            @testset "stationarity against exact enumeration (stacked triangular)" begin
                # calibration-begin
                # Offset (2, 2, 3) cell, 24 sites, 5 particles (42504
                # configurations), in-plane nearest-neighbour coupling in
                # shell 1 and the adjacent-layer coupling in shell 2. A seeded
                # fixed-N nested-sampling ladder with mixed local and cluster
                # moves must reproduce the exact canonical mean energy at two
                # temperatures; a local-swap-only run at the same step budget
                # is the control at the same gate.
                st_N = 5
                st_lattice() = offset_tri(2, 2, 3; comps=[collect(1:st_N)])
                st_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.006], u"eV")
                st_kB = 8.617333262e-5
                st_betas = [1 / (st_kB * 150.0), 1 / (st_kB * 400.0)]
                st_K = 200
                st_nsteps = Int64(3000)
                function st_run(seed, clusters_freq)
                    Random.seed!(seed)
                    walkers = LatticeWalker{1}[]
                    for _ in 1:st_K
                        lat = st_lattice()
                        lat.components[1] .= false
                        lat.components[1][randperm(num_sites(lat))[1:st_N]] .= true
                        push!(walkers, LatticeWalker(lat, energy=0.0u"eV", iter=0))
                    end
                    ls = LatticeGasWalkers(walkers, st_ham)
                    params = LatticeNestedSamplingParameters(mc_steps=40,
                                                             allowed_fail_count=10^9)
                    routine = MCMixedMoves(walks_freq=1, clusters_freq=clusters_freq,
                                           initial_cluster_p=0.3,
                                           cluster_adjust_interval=50)
                    tag = "t_st_$(seed)_$(clusters_freq)"
                    save = SaveEveryN("$(tag).csv", "$(tag).traj", "$(tag).ls",
                                      10^6, 10^6, 10^6)
                    df, _, params_out = nested_sampling(ls, params, st_nsteps,
                                                        routine, save)
                    rm.(["$(tag).csv", "$(tag).traj", "$(tag).ls"], force=true)
                    w = ωᵢ(df.iter, st_K)
                    return [internal_energy(β, w, df.emax) for β in st_betas], params_out
                end
                st_df_ex, _ = exact_enumeration(st_lattice(), st_ham)
                st_E_ex = [e.val for e in st_df_ex.energy]
                st_U_ex = [internal_energy(β, ones(length(st_E_ex)), st_E_ex)
                           for β in st_betas]
                # calibration-end
                st_U_mix, st_p_mix = st_run(7301, 1)
                st_U_ctl, _ = st_run(7301, 0)
                # Gates: 3x the maximum three-seed deviation per temperature
                # over both routines (seeds 7301, 7302, 7303; the shipped seed
                # is the first). Calibration |U_NS - U_exact| in eV:
                #   mixed   150 K: 2.1858e-3, 9.502e-4, 8.660e-4
                #           400 K: 6.152e-4, 1.239e-4, 1.4213e-3
                #   control 150 K: 5.899e-4, 3.345e-4, 9.251e-4
                #           400 K: 4.443e-4, 1.093e-4, 1.0436e-3
                # Maxima 2.1858e-3 and 1.4213e-3, both from the mixed run;
                # U_exact = -0.26379 eV and -0.24833 eV.
                st_gate = [0.0066, 0.0043]
                for t in 1:2
                    @test abs(st_U_mix[t] - st_U_ex[t]) < st_gate[t]
                    @test abs(st_U_ctl[t] - st_U_ex[t]) < st_gate[t]
                end
                @test any(>(0.0), st_p_mix.cluster_accept_history)
            end

            @testset "drivers on a stacked offset cell" begin
                drv_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.006], u"eV")
                drv_cleanup(tag) = rm.(["$(tag).csv", "$(tag).traj", "$(tag).ls"],
                                       force=true)
                drv_save(tag) = SaveEveryN("$(tag).csv", "$(tag).traj", "$(tag).ls",
                                           10^6, 10^6, 10^6)

                # Fixed-N nested sampling with mixed local and cluster moves
                Random.seed!(7401)
                walkers = LatticeWalker{1}[]
                for _ in 1:20
                    lat = offset_tri(4, 4, 3; comps=[Int[]])
                    lat.components[1][randperm(num_sites(lat))[1:24]] .= true
                    push!(walkers, LatticeWalker(lat, energy=0.0u"eV", iter=0))
                end
                ls = LatticeGasWalkers(walkers, drv_ham)
                params = LatticeNestedSamplingParameters(mc_steps=40,
                                                         allowed_fail_count=10^9)
                routine = MCMixedMoves(walks_freq=1, clusters_freq=1,
                                       cluster_adjust_interval=10)
                df, ls_out, p_out = nested_sampling(ls, params, Int64(300), routine,
                                                    drv_save("t_drv_fn"))
                drv_cleanup("t_drv_fn")
                @test df isa DataFrame
                @test nrow(df) > 0
                @test length(ls_out.walkers) == 20
                @test all(sum(w.configuration.components[1]) == 24
                          for w in ls_out.walkers)
                @test any(>(0.0), p_out.cluster_accept_history)

                # Ideal-gas-referenced grand-canonical driver with cluster moves
                Random.seed!(7402)
                gw = [LatticeWalker(offset_tri(4, 4, 3; comps=[Int[]]),
                                    energy=0.0u"eV", iter=0) for _ in 1:20]
                gls = LatticeGasWalkers(gw, drv_ham; assign_energy=false)
                gparams = IdealGasReferencedGCNSParameters(mc_steps=40,
                                                           reference_fugacity=0.5)
                groutine = MCGrandCanonicalMoves(p_move=0.5, p_insert=0.25,
                                                 clusters_freq=1, swaps_freq=1)
                gdf, _, gp = ideal_gas_referenced_nested_sampling(
                    gls, gparams, Int64(300), groutine, drv_save("t_drv_ig"))
                drv_cleanup("t_drv_ig")
                @test gdf isa DataFrame
                @test nrow(gdf) > 0
                @test gp.move_stats[:cluster_attempted] > 0
                @test gp.move_stats[:cluster_accepted] > 0
            end
        end

    end

    @testset "copy-free proposal and revert" begin
        using Random
        using Unitful

        cf_lattice() = MLattice{1,SquareLattice}(lattice_constant=1.0,
            basis=[(0.0, 0.0, 0.0)], supercell_dimensions=(4, 4, 1),
            periodicity=(true, true, false), cutoff_radii=[1.1],
            components=[[false for _ in 1:16]], adsorptions=:full)
        cf_ham() = GenericLatticeHamiltonian(-0.04, [-0.01], u"eV")
        cf_tri() = SLattice{TriangularLattice}(
            supercell_dimensions=(6, 4, 1),
            cutoff_radii=[1.1],
            components=[[1, 5, 10, 20, 30, 40]],
            adsorptions=:full)

        @testset "configuration identity through every kernel" begin
            Random.seed!(91001)
            lat = cf_lattice()
            for i in 1:16
                lat.components[1][i] = rand() < 0.5
            end
            wk = LatticeWalker(lat, energy=interacting_energy(lat, cf_ham()),
                               iter=0)
            cfg0 = wk.configuration
            MC_random_walk!(50, wk, cf_ham(), 1.0e3; energy_perturb=1e-9)
            @test wk.configuration === cfg0
            MC_cluster_walk!(20, wk, cf_ham(), 1.0e3, 0.3; energy_perturb=1e-9)
            @test wk.configuration === cfg0
            MC_grand_canonical_walk!(100, wk, cf_ham(), 1.0e3, 0.0;
                p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
                clusters_freq=2, swaps_freq=2, cluster_p=0.3)
            @test wk.configuration === cfg0
        end

        @testset "forced-reject revert exactness" begin
            # A ceiling below every reachable energy rejects every proposal;
            # occupancy and the stored energy must be exactly unchanged.
            Random.seed!(91002)
            lat = cf_lattice()
            for i in 1:16
                lat.components[1][i] = rand() < 0.5
            end
            e0 = interacting_energy(lat, cf_ham())
            wk = LatticeWalker(lat, energy=e0, iter=0)
            occ0 = copy(lat.components[1])
            low = -1.0e3

            a1, r1, _ = MC_random_walk!(60, wk, cf_ham(), low;
                                        energy_perturb=1e-9)
            @test a1 == false && r1 == 0.0
            @test wk.configuration.components[1] == occ0
            @test wk.energy == e0

            a2, r2, _ = MC_cluster_walk!(30, wk, cf_ham(), low, 0.3;
                                         energy_perturb=1e-9)
            @test a2 == false && r2 == 0.0
            @test wk.configuration.components[1] == occ0
            @test wk.energy == e0

            a3, r3, _, _, _, ms = MC_grand_canonical_walk!(200, wk, cf_ham(),
                low, 0.0;
                p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
                clusters_freq=2, swaps_freq=2, cluster_p=0.3)
            @test a3 == false && r3 == 0.0
            @test wk.configuration.components[1] == occ0
            @test wk.energy == e0
            # every default channel proposed under the forced-reject ceiling
            @test ms.swap_attempted > 0
            @test ms.cluster_attempted > 0
            @test ms.insert_uniform_attempted > 0
            @test ms.delete_attempted > 0

            # composite channel: biased insertion and the inlined deletion
            # revert through the same path
            a4, _, _, _, _, ms4 = MC_grand_canonical_walk!(200, wk, cf_ham(),
                low, 0.0;
                p_move=0.4, p_insert=0.3, z0=1.0, energy_perturb=1e-9,
                swaps_freq=1, p_bias=0.4)
            @test a4 == false
            @test wk.configuration.components[1] == occ0
            @test wk.energy == e0
            @test ms4.insert_biased_attempted > 0

            # triangular cluster revert exercises the recorded-pair replay
            # on the second geometry
            Random.seed!(91003)
            tri = cf_tri()
            et = interacting_energy(tri, cf_ham())
            wt = LatticeWalker(tri, energy=et, iter=0)
            occt = copy(tri.components[1])
            at, rt, _ = MC_cluster_walk!(30, wt, cf_ham(), low, 0.3;
                                         energy_perturb=1e-9)
            @test at == false && rt == 0.0
            @test wt.configuration.components[1] == occt
            @test wt.energy == et
        end

        @testset "multi-component forced-reject revert" begin
            # The C > 1 revert path replays through the empty/occupied and
            # cross-component exchanges (both involutions): a forced-reject
            # walk leaves every component vector exactly unchanged and the
            # configuration object identity intact
            Random.seed!(91008)
            ml = MLattice{2,SquareLattice}(components=[[1, 3, 6], [2, 4, 7]])
            mlham = MLatticeHamiltonian(2,
                [cf_ham(), GenericLatticeHamiltonian(-0.02, [-0.005], u"eV"),
                 cf_ham()])
            wm = LatticeWalker(ml, energy=interacting_energy(ml, mlham),
                               iter=0)
            comps0 = [copy(v) for v in ml.components]
            e0m = wm.energy
            cfg0m = wm.configuration
            am, rm, _ = MC_random_walk!(80, wm, mlham, -1.0e3;
                                        energy_perturb=1e-9)
            @test am == false && rm == 0.0
            @test wm.configuration === cfg0m
            @test wm.configuration.components[1] == comps0[1]
            @test wm.configuration.components[2] == comps0[2]
            @test wm.energy == e0m
        end

        @testset "cluster record keyword" begin
            # Recorded pairs re-applied restore the original configuration
            Random.seed!(91004)
            lat = cf_lattice()
            for i in 1:16
                lat.components[1][i] = rand() < 0.5
            end
            ref = deepcopy(lat.components[1])
            pairs = Tuple{Int,Int}[]
            geometric_cluster_swap!(lat, 0.4; record=pairs)
            MonteCarloMoves._apply_cluster_pairs!(lat, pairs)
            @test lat.components[1] == ref

            Random.seed!(91007)
            tri = cf_tri()
            reft = deepcopy(tri.components[1])
            pt = Tuple{Int,Int}[]
            geometric_cluster_swap!(tri, 0.4; record=pt)
            MonteCarloMoves._apply_cluster_pairs!(tri, pt)
            @test tri.components[1] == reft

            # The nothing default changes neither the configuration nor the
            # random stream: same-seed A/B with and without a record vector
            Random.seed!(91005)
            latA = cf_lattice()
            for i in 1:16
                latA.components[1][i] = rand() < 0.5
            end
            latB = deepcopy(latA)
            Random.seed!(91006)
            geometric_cluster_swap!(latA, 0.4)
            nextA = rand()
            Random.seed!(91006)
            geometric_cluster_swap!(latB, 0.4; record=Tuple{Int,Int}[])
            nextB = rand()
            @test latA.components[1] == latB.components[1]
            @test nextA == nextB

            # The generic-geometry guard keeps its descriptive error with the
            # keyword passed through
            gen = MLattice{1,GenericLattice}(
                [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
                [(0.0, 0.0, 0.0)],
                (3, 1, 1),
                (false, false, false),
                [1.1],
                [[true, false, false]],
                fill(false, 3))
            @test_throws ArgumentError geometric_cluster_swap!(gen, 0.3;
                record=Tuple{Int,Int}[])
        end
    end
end
