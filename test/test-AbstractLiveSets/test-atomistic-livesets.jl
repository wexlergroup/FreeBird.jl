@testset "Atomistic live sets tests" begin

    @testset "LJAtomWalkers & GuptaAtomWalkers struct and functions tests" begin

        @testset "Constructor with energy assignment" begin

            # Create 10 random configurations for testing
            configs = generate_initial_configs(10, 10.0, 10)
            # Create walkers from configurations
            walkers = AtomWalker.(configs)
            # Create a Lennard-Jones potential
            lj = LJParameters()
            # Create a Gupta potential
            gp = GuptaParameters()
            # Create live sets to test
            lj_walkers = LJAtomWalkers(walkers, lj)
            gp_walkers = GuptaAtomWalkers(walkers, gp)

            # Test type and structure
            @test lj_walkers isa LJAtomWalkers
            @test gp_walkers isa GuptaAtomWalkers

            for ls in [lj_walkers, gp_walkers] 
                @test length(ls.walkers) == 10
                @test ls.potential !== nothing
                
                # Test energy assignment
                @test all(w.energy != 0.0u"eV" for w in ls.walkers)
                @test all(w.energy_frozen_part == 0.0u"eV" for w in ls.walkers)
            end

        end

        @testset "Constructor without energy assignment" begin

            # Create 10 random configurations for testing
            configs = generate_initial_configs(10, 10.0, 10)
            # Create walkers from configurations
            walkers = AtomWalker.(configs)
            # Create a Lennard-Jones potential
            lj = LJParameters()
            # Create a Gupta potential
            gp = GuptaParameters()

            # Create live sets to test
            lj_walkers = LJAtomWalkers(walkers, lj, assign_energy=false)
            gp_walkers = GuptaAtomWalkers(walkers, gp, assign_energy=false)

            # Verify energies weren't assigned
            for ls in [lj_walkers, gp_walkers] 
                @test all(w.energy == 0.0u"eV" for w in ls.walkers)
                @test all(w.energy_frozen_part == 0.0u"eV" for w in ls.walkers)
            end
        end
    end

    @testset "LJSurfaceWalkers struct and functions tests" begin
        
        # Create 10 random walkers
        walkers = AtomWalker.(generate_initial_configs(10, 10.0, 10))

        # Create a fictitious surface with 4 H atoms at z=0 plane
        surf_list = [:H => [0.0, 0.0, 0.0], :H => [0.0, 0.5, 0.0], :H => [0.5, 0.0, 0.0], :H => [0.5, 0.5, 0.0]]
        box = walkers[1].configuration.cell.cell_vectors # extract box from one of the walkers

        # Create AtomWalker for the (frozen) surface
        surface = AtomWalker(FastSystem(periodic_system(surf_list, box, fractional=true)); freeze_species=[:H])

        # Composite LJ parameters, should have two components: one for free particles and one for surface
        ljs = CompositeParameterSets(2, [LJParameters(epsilon=i) for i in 1:3]) # 2 components, 3 parameter sets

        # Pre-calculate frozen part energy for the surface, not that it matters since it's frozen, making it non-zero here for testing purposes
        surface.energy_frozen_part = interacting_energy(surface.configuration, ljs.param_sets[2,2])

        @testset "Constructor with energy assignment" begin

            # Create a live set, single-threaded/single-process by default
            lj_walkers = LJSurfaceWalkers(walkers, ljs, surface; assign_energy=true)

            # Test type and structure
            @test lj_walkers isa LJSurfaceWalkers
            @test length(lj_walkers.walkers) == 10
            @test lj_walkers.potential === ljs
            @test lj_walkers.surface === surface
            
            # Test energy assignment
            @test all(w.energy != 0.0u"eV" for w in lj_walkers.walkers)
            @test all(w.energy_frozen_part == interacting_energy(surface.configuration, ljs.param_sets[2,2]) for w in lj_walkers.walkers)

            # Test thread safety
            threaded_lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, :threads)
            @test threaded_lj_walkers isa LJSurfaceWalkers
            @test length(threaded_lj_walkers.walkers) == 10
            @test threaded_lj_walkers.potential === ljs
            @test threaded_lj_walkers.surface === surface
            @test all(threaded_lj_walkers.walkers[i].energy == lj_walkers.walkers[i].energy for i in 1:10)
            @test all(threaded_lj_walkers.walkers[i].energy_frozen_part == lj_walkers.walkers[i].energy_frozen_part for i in 1:10)

            # Test distributed safety
            distributed_lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, :distributed)
            @test distributed_lj_walkers isa LJSurfaceWalkers
            @test length(distributed_lj_walkers.walkers) == 10
            @test distributed_lj_walkers.potential === ljs
            @test distributed_lj_walkers.surface === surface
            @test all(distributed_lj_walkers.walkers[i].energy == lj_walkers.walkers[i].energy for i in 1:10)
            @test all(distributed_lj_walkers.walkers[i].energy_frozen_part == lj_walkers.walkers[i].energy_frozen_part for i in 1:10)
        end

        @testset "Constructor without energy assignment" begin
            # Reset walker energies
            walkers = AtomWalker.(generate_initial_configs(10, 10.0, 10))

            lj_walkers = LJSurfaceWalkers(walkers, ljs, surface, assign_energy=false)

            # Verify energies weren't assigned
            @test all(w.energy == 0.0u"eV" for w in lj_walkers.walkers)
            @test all(w.energy_frozen_part == 0.0u"eV" for w in lj_walkers.walkers)
        end
    end

end