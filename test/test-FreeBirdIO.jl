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

        save_free_part_every_n = SaveFreePartEveryN()
        @test save_every_n.df_filename == "output_df.csv"
        @test save_every_n.wk_filename == "output.traj.extxyz"
        @test save_every_n.ls_filename == "output.ls.extxyz"
        @test save_every_n.n_traj == 100
        @test save_every_n.n_snap == 1000
    end

    @testset "extract free part" begin
        sys = FreeBirdIO.generate_multi_type_random_starting_config(10.0, [1, 2, 3], particle_types=[:H, :O, :H])
        walker = AtomWalker(sys; freeze_species=[:H])
        free = FreeBirdIO.extract_free_par(walker)
        @test length(free.configuration) == 2
        @test all([!frozen for frozen in free.frozen])
        @test all(ChemicalSpecies(:O) in free.configuration.species)
    end


end