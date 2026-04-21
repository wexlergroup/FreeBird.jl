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

    end
end
