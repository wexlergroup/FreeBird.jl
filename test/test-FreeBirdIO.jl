@testset "FreeBirdIO Tests" begin
    
    @testset "set_pbc" begin
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
        boundary_conditions = (true, true, true)
        sys = FlexibleSystem([Atom(:H, [0, 0, 1.]u"Å")], box, boundary_conditions)
        at = Atoms(sys)
        pbc = [true, false, true]
        flex_sys = FreeBirdIO.set_pbc(at, pbc)
        @test periodicity(flex_sys) == Tuple(pbc)
    end

    @testset "convert_system_to_walker" begin
        at = FlexibleSystem(FreeBirdIO.generate_random_starting_config(100.0, 1))
        walker = convert_system_to_walker(at, false)
        @test walker.energy == 0.0u"eV"
        @test walker.iter == 0
        @test walker.list_num_par == [1]
        @test walker.frozen == [false]
        @test walker.energy_frozen_part == 0.0u"eV"
    end

    @testset "convert_walker_to_system" begin
        walker = AtomWalker(FreeBirdIO.generate_random_starting_config(100.0, 1))
        sys = convert_walker_to_system(walker)
        @test sys isa AbstractSystem
    end

    @testset "append_system" begin
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
        boundary_conditions = (true, true, true)
        at1 = FlexibleSystem([Atom(:H, [0, 0, 1.]u"Å")], box, boundary_conditions)
        at2 = FlexibleSystem([Atom(:H, [0, 0, 3.]u"Å")], box, boundary_conditions)
        new_sys = append_system(at1, at2)
        @test new_sys isa FastSystem
        @test length(new_sys.position) == 2
    end

    @testset "generate_random_starting_config" begin
        volume_per_particle = 10.0
        num_particle = 5
        sys = FreeBirdIO.generate_random_starting_config(volume_per_particle, num_particle)
        @test sys isa FastSystem
        @test length(sys.position) == num_particle
    end

    @testset "generate_multi_type_random_starting_config" begin
        volume_per_particle = 10.0
        num_particle = [3, 2]
        particle_types = [:H, :O]
        sys = FreeBirdIO.generate_multi_type_random_starting_config(volume_per_particle, num_particle; particle_types=particle_types)
        @test sys isa FastSystem
        @test length(sys.position) == sum(num_particle)
    end

    @testset "generate_initial_configs" begin
        num_walkers = 3
        volume_per_particle = 10.0
        num_particle = 5
        configs = generate_initial_configs(num_walkers, volume_per_particle, num_particle)
        @test length(configs) == num_walkers
        @test all(config -> config isa FastSystem, configs)
    end

    @testset "DataSavingStrategy structs" begin
        save_every_n = SaveEveryN()
        @test save_every_n.df_filename == "output_df.csv"
        @test save_every_n.wk_filename == "output.traj.extxyz"
        @test save_every_n.ls_filename == "output.ls.extxyz"
        @test save_every_n.n_traj == 100
        @test save_every_n.n_snap == 1000

        save_every_n_custom = SaveEveryN("custom_df.csv", "custom.traj.extxyz", "custom.ls.extxyz", 200, 2000, 200)

        @test save_every_n_custom.df_filename == "custom_df.csv"
        @test save_every_n_custom.wk_filename == "custom.traj.extxyz"
        @test save_every_n_custom.ls_filename == "custom.ls.extxyz"
        @test save_every_n_custom.n_traj == 200
        @test save_every_n_custom.n_snap == 2000
        @test save_every_n_custom.n_info == 200

        save_free_part_every_n = SaveFreePartEveryN()
        @test save_every_n.df_filename == "output_df.csv"
        @test save_every_n.wk_filename == "output.traj.extxyz"
        @test save_every_n.ls_filename == "output.ls.extxyz"
        @test save_every_n.n_traj == 100
        @test save_every_n.n_snap == 1000
        @test save_every_n.n_info == 1

        save_free_part_every_n_custom = SaveFreePartEveryN("custom_df.csv", "custom.traj.extxyz", "custom.ls.extxyz", 200, 2000, 200)
        @test save_every_n_custom.df_filename == "custom_df.csv"
        @test save_every_n_custom.wk_filename == "custom.traj.extxyz"
        @test save_every_n_custom.ls_filename == "custom.ls.extxyz"
        @test save_every_n_custom.n_traj == 200
        @test save_every_n_custom.n_snap == 2000
    end

    @testset "extract free part" begin
        sys = FreeBirdIO.generate_multi_type_random_starting_config(10.0, [1, 2, 3], particle_types=[:H, :O, :H])
        walker = AtomWalker(sys; freeze_species=[:H])
        free = FreeBirdIO.extract_free_par(walker)
        @test length(free.configuration) == 2
        @test all([!frozen for frozen in free.frozen])
        @test all(ChemicalSpecies(:O) in free.configuration.species)
    end

    @testset "save and read configurations" begin
        volume_per_particle = 10.0
        num_particle = [3, 2]
        particle_types = [:H, :O]
        sys = FreeBirdIO.generate_multi_type_random_starting_config(volume_per_particle, num_particle; particle_types=particle_types)
        walker = AtomWalker(sys)
        FreeBirdIO.write_single_walker("test.traj.extxyz", walker)
        @test isfile("test.traj.extxyz")
        readin_config = FreeBirdIO.read_single_config("test.traj.extxyz")
        @test readin_config isa FlexibleSystem
        @test length(readin_config) == sum(num_particle)
        @test periodicity(readin_config) == (true, true, true) # which is wrong
        readin_config = FreeBirdIO.read_single_config("test.traj.extxyz", pbc="FFF")
        @test periodicity(readin_config) == (false, false, false) # which is correct
        rm("test.traj.extxyz", force=true)

        walkers  = [AtomWalker(FreeBirdIO.generate_multi_type_random_starting_config(volume_per_particle, num_particle; particle_types=particle_types)) for _ in 1:3]
        write_walkers("walkers.traj.extxyz", walkers)
        @test isfile("walkers.traj.extxyz")
        readin_walkers = read_walkers("walkers.traj.extxyz"; resume=false)
        @test length(readin_walkers) == 3
        @test all(walker -> walker isa AtomWalker, readin_walkers)

        readin_single = read_single_walker("walkers.traj.extxyz"; resume=false)
        @test readin_single isa AtomWalker
        @test periodicity(readin_single.configuration) == (true, true, true) # which is wrong
        readin_single = read_single_walker("walkers.traj.extxyz", pbc="FFF"; resume=false)
        @test periodicity(readin_single.configuration) == (false, false, false) # which is correct

        rm("walkers.traj.extxyz", force=true)
        
    end
        


end
@testset "Overlap-aware periodic starting configurations" begin
    # Minimum-image pairwise distance under the generated system's own metric
    function _test_min_dist(sys)
        L = ustrip(u"Å", cell_vectors(sys)[1][1])
        pbc = periodicity(sys)
        n = length(sys)
        dmin = Inf
        for i in 1:(n - 1), j in (i + 1):n
            d2 = 0.0
            for k in 1:3
                δ = ustrip(u"Å", position(sys, i)[k] - position(sys, j)[k])
                if pbc[k]
                    δ -= L * round(δ / L)
                end
                d2 += δ^2
            end
            dmin = min(dmin, sqrt(d2))
        end
        return dmin
    end

    @testset "Default path is stream-neutral" begin
        Random.seed!(31415)
        sys_new = FreeBirdIO.generate_random_starting_config(100.0, 5)
        Random.seed!(31415)
        box_length = (100.0 * 5)^(1 / 3)
        legacy_positions = [[rand(), rand(), rand()] .* box_length * u"Å" for _ in 1:5]
        @test all(position(sys_new, i) == SVector{3}(legacy_positions[i]) for i in 1:5)
        @test periodicity(sys_new) == (false, false, false)
    end

    @testset "Periodicity propagation" begin
        sys = FreeBirdIO.generate_random_starting_config(100.0, 3; periodicity=(true, true, false))
        @test periodicity(sys) == (true, true, false)
        sys_p = FreeBirdIO.generate_random_starting_config(100.0, 3; periodicity=(true, true, true))
        @test periodicity(sys_p) == (true, true, true)
    end

    @testset "Separation invariant at a dense point" begin
        # 24 particles at 15.625 A^3 each (box 7.211 A): packing fraction 0.257
        # at the 2.125 A separation, safely under the sequential-insertion regime
        Random.seed!(2718)
        sys = FreeBirdIO.generate_random_starting_config(15.625, 24;
            periodicity=(true, true, true), min_separation=2.125, max_attempts=5000)
        @test length(sys) == 24
        @test _test_min_dist(sys) >= 2.125
        Random.seed!(2718)
        sys_m = FreeBirdIO.generate_random_starting_config(15.625, 24;
            periodicity=(true, true, false), min_separation=2.125, max_attempts=5000)
        @test _test_min_dist(sys_m) >= 2.125
    end

    @testset "Infeasible packing throws" begin
        Random.seed!(1618)
        @test_throws ArgumentError FreeBirdIO.generate_random_starting_config(15.625, 24;
            periodicity=(true, true, true), min_separation=4.0, max_attempts=200)
    end

    @testset "Lattice starting configuration" begin
        sys1 = FreeBirdIO.generate_lattice_starting_config(12.5, 100)
        sys2 = FreeBirdIO.generate_lattice_starting_config(12.5, 100)
        @test length(sys1) == 100
        @test periodicity(sys1) == (true, true, true)
        # deterministic at jitter = 0
        @test all(position(sys1, i) == position(sys2, i) for i in 1:100)
        # separation bound: at least the grid spacing (100 sites on a 5^3 grid)
        @test _test_min_dist(sys1) >= 12.5 / 5 - 1e-9
        # partial fill keeps the bound
        sys3 = FreeBirdIO.generate_lattice_starting_config(12.5, 30)
        @test length(sys3) == 30
        @test _test_min_dist(sys3) >= 12.5 / ceil(Int, cbrt(30)) - 1e-9
        # jitter: seeded reproducibility and the documented separation bound
        Random.seed!(999)
        sysj1 = FreeBirdIO.generate_lattice_starting_config(12.5, 100; jitter=0.2)
        Random.seed!(999)
        sysj2 = FreeBirdIO.generate_lattice_starting_config(12.5, 100; jitter=0.2)
        @test all(position(sysj1, i) == position(sysj2, i) for i in 1:100)
        @test _test_min_dist(sysj1) >= 12.5 / 5 - 2 * sqrt(3) * 0.2 - 1e-9
        @test_throws ArgumentError FreeBirdIO.generate_lattice_starting_config(12.5, 0)
    end
end
