@testset "nested_sampling.jl tests" begin
    @testset "NestedSamplingParameters struct tests" begin
        params = NestedSamplingParameters(
            1000,    # mc_steps
            0.1,     # initial_step_size
            0.1,     # step_size
            0.01,    # step_size_lo
            1.0,     # step_size_up
            (0.5,0.5),  # acceptance ratio range
            0,       # fail_count
            100,     # allowed_fail_count
            1e-12,   # energy_perturbation
            1234     # random_seed
        )
        
        @test params isa SamplingSchemes.SamplingParameters
        @test params.mc_steps == 1000
        @test params.initial_step_size == 0.1
        @test params.step_size == 0.1
        @test params.step_size_lo == 0.01
        @test params.step_size_up == 1.0
        @test params.accept_range == (0.5,0.5)
        @test params.fail_count == 0
        @test params.allowed_fail_count == 100
        @test params.energy_perturbation == 1e-12
        @test params.random_seed == 1234
        
        # Test mutability
        params.step_size = 0.2
        params.fail_count = 1
        @test params.step_size == 0.2
        @test params.fail_count == 1

        params = NestedSamplingParameters() # Default values

        @test params.mc_steps == 200
        @test params.initial_step_size == 0.01
        @test params.step_size == 0.1
        @test params.step_size_lo == 1e-6
        @test params.step_size_up == 1.0
        @test params.accept_range == (0.25,0.75)
        @test params.fail_count == 0
        @test params.allowed_fail_count == 100
        @test params.energy_perturbation == 1e-12
        @test params.random_seed == 1234

    end
    

    @testset "LatticeNestedSamplingParameters struct tests" begin
        params = LatticeNestedSamplingParameters(
            mc_steps=1000,
            energy_perturbation=0.1,
            fail_count=0,
            allowed_fail_count=100,
            random_seed=1234
        )
        
        @test params isa SamplingSchemes.SamplingParameters
        @test params.mc_steps == 1000
        @test params.energy_perturbation == 0.1
        @test params.fail_count == 0
        @test params.allowed_fail_count == 100
        @test params.random_seed == 1234
        
        # Test mutability
        params.fail_count = 1
        @test params.fail_count == 1

        params = LatticeNestedSamplingParameters() # Default values

        @test params.mc_steps == 100
        @test params.energy_perturbation == 1e-12
        @test params.fail_count == 0
        @test params.allowed_fail_count == 10
        @test params.random_seed == 1234

    end
    

    @testset "MCRoutine types tests" begin
        @test MCRandomWalkMaxE() isa MCRoutine
        @test MCRandomWalkClone() isa MCRoutine
        @test MCNewSample() isa MCRoutine
        @test MCMixedMoves(5,1) isa MCRoutine
        @test MCRandomWalkClone() isa MCRoutine
        @test MCRandomWalkMaxEParallel() isa MCRoutine
    end


    @testset "AtomWalkers functions tests" begin
        # Create periodic box and LJ parameters
        box = [[10.0u"Å", 0u"Å", 0u"Å"],
                [0u"Å", 10.0u"Å", 0u"Å"],
                [0u"Å", 0u"Å", 10.0u"Å"]]

        lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

        # Test with 4 atoms: 2 H and 2 O atoms
        coor_list = [:H => [0.2, 0.5, 0.5],
                    :H => [0.4, 0.5, 0.5],
                    :O => [0.6, 0.5, 0.5],
                    :O => [0.8, 0.5, 0.5]]
        at = FastSystem(periodic_system(coor_list, box, fractional=true))

        # Create multiple walkers with different configurations
        walkers = [
            AtomWalker(at),  # No frozen atoms
            AtomWalker(at, freeze_species=[:O]),  # O atoms frozen
            AtomWalker(at, freeze_species=[:H, :O])  # All frozen
        ]

        # Create LJAtomWalkers instance
        lj_set = LJAtomWalkers(walkers, lj)

        @testset "sort_by_energy function tests" begin
            # Store original energies
            original_energies = [w.energy for w in lj_set.walkers]
            
            # Sort walkers
            sorted_set = sort_by_energy!(lj_set)
            
            # Test sorting is in descending order
            sorted_energies = [w.energy for w in sorted_set.walkers]
            @test issorted(sorted_energies, rev=true)
            
            # Test that all original energies are present
            @test sort(original_energies) == sort(sorted_energies)
            
            # Test that the original set was modified (not a copy)
            @test lj_set === sorted_set
            
            # Test sorting is stable for equal energies
            equal_energy_walkers = [
                AtomWalker(at),
                AtomWalker(at)
            ]
            equal_energy_set = LJAtomWalkers(equal_energy_walkers, lj)
            original_order = [w.configuration for w in equal_energy_set.walkers]
            sort_by_energy!(equal_energy_set)
            final_order = [w.configuration for w in equal_energy_set.walkers]
            @test original_order == final_order
        end

        @testset "update_iter function tests" begin
            # Check initial iterations
            @test all(w.iter == 0 for w in lj_set.walkers)
            
            # Update iterations once
            SamplingSchemes.update_iter!(lj_set)
            @test all(w.iter == 1 for w in lj_set.walkers)
            
            # Update iterations multiple times
            for i in 2:5
                SamplingSchemes.update_iter!(lj_set)
                @test all(w.iter == i for w in lj_set.walkers)
            end
            
            # Test with empty walker set
            empty_set = LJAtomWalkers(AtomWalker{typeof(at)}[], lj)
            SamplingSchemes.update_iter!(empty_set)
            @test isempty(empty_set.walkers)
        end
    end


    @testset "nested_sampling_step functions tests" begin
        @testset "AtomWalkers nested_sampling_step!" begin
            # Set up
            box = [[10.0u"Å", 0u"Å", 0u"Å"],
                    [0u"Å", 10.0u"Å", 0u"Å"],
                    [0u"Å", 0u"Å", 10.0u"Å"]]

            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
            ljs = CompositeLJParameters(3, [lj for _ in 1:6])
            
            coor_list = [:H => [0.2, 0.5, 0.5],
                        :H => [0.4, 0.5, 0.5],
                        :O => [0.6, 0.5, 0.5]]
            at = FastSystem(periodic_system(coor_list, box, fractional=true))

            surf_list = [:H => [0.0, 0.0, 0.0],
                     :H => [0.0, 0.5, 0.0],
                     :H => [0.5, 0.0, 0.0],
                     :H => [0.5, 0.5, 0.0]]

            surface = AtomWalker(FastSystem(periodic_system(surf_list, box, fractional=true)); freeze_species=[:H])
            surface.energy_frozen_part = interacting_energy(surface.configuration, lj)
            
            walkers = [AtomWalker(at) for _ in 1:3]
            liveset_at = LJAtomWalkers(walkers, lj)
            liveset_surf = LJSurfaceWalkers(walkers, ljs, surface; assign_energy=true)
            
            ns_params = NestedSamplingParameters(
                1000,    # mc_steps
                0.1,     # initial_step_size
                0.1,     # step_size
                0.01,    # step_size_lo
                1.0,     # step_size_up
                (0.5,0.5),  # acceptance ratio range
                0,       # fail_count
                100,     # allowed_fail_count
                1e-12,   # energy_perturbation
                1234     # random_seed
            )
            
            for liveset in [liveset_at, liveset_surf]
                @testset "MCRandomWalkMaxE" begin
                    mc_routine = MCRandomWalkMaxE()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,typeof(0.0u"eV")}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end

                @testset "MCRandomWalkMaxEParallel" begin
                    mc_routine = MCRandomWalkMaxEParallel()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,typeof(0.0u"eV")}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end
        
                @testset "MCRandomWalkClone" begin
                    mc_routine = MCRandomWalkClone()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,typeof(0.0u"eV")}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end


                @testset "MCRandomWalkCloneParallel" begin
                    mc_routine = MCRandomWalkCloneParallel()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,typeof(0.0u"eV")}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end



                @testset "MCRandomWalkClone 2D" begin
                    mc_routine = MCRandomWalkClone(dims=[1,2])
                    @test mc_routine.dims == [1,2]
                    @test length(mc_routine.dims) == 2

                    original_ls = deepcopy(liveset)

                    if !(liveset isa LJSurfaceWalkers) # TODO: LJSurfaceWalkers does not support 2D walks yet
                        iter, emax, updated_liveset, updated_params = nested_sampling_step!(deepcopy(original_ls), ns_params, mc_routine)
                        
                        @test iter isa Union{Missing,Int}
                        @test emax isa Union{Missing,typeof(0.0u"eV")}
                        @test length(updated_liveset.walkers) == length(original_ls.walkers)
                        @test updated_params.fail_count >= 0
                    end

                end

                @testset "MCRandomWalkClone 1D" begin
                    mc_routine = MCRandomWalkClone(dims=[1])
                    @test mc_routine.dims == [1]
                    @test length(mc_routine.dims) == 1

                    @test_throws ErrorException begin
                        nested_sampling_step!(liveset, ns_params, mc_routine)
                    end

                end

                @testset "MCMixedMoves" begin
                    mc_routine = MCMixedMoves(5, 1)

                    if !(liveset isa LJSurfaceWalkers) # TODO: LJSurfaceWalkers does not support mixed moves yet
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                        @test iter isa Union{Missing,Int}
                        @test emax isa Union{Missing,typeof(0.0u"eV")}
                        @test length(updated_liveset.walkers) == length(liveset.walkers)
                        @test updated_params.fail_count >= 0
                    end
                end

                @testset "Unsupported MCRoutine" begin
                    mc_routine = MCNewSample()
                    struct Unsupported <: MCRoutine end
                    @test_throws ErrorException begin
                        nested_sampling_step!(liveset, ns_params, Unsupported())
                    end
                end
            end
        end

        @testset "LatticeGasWalkers nested_sampling_step!" begin
            # Setup for LatticeGasWalkers
            square_lattice = MLattice{1,SquareLattice}(
                    lattice_constant=1.0,                   # Unit cell size
                    basis=[(0.0, 0.0, 0.0)],                # Single atom per unit cell
                    supercell_dimensions=(4, 4, 1),         # 4x4x1 supercell
                    periodicity=(true, true, false),        # Periodic in x,y but not z
                    cutoff_radii=[1.1, 1.5],                # Neighbor cutoff distances
                    components=:equal,                      # Split into equal components
                    adsorptions=:full                       # All sites are adsorption sites
                )

            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

            s_walker = LatticeWalker(
                square_lattice,
                energy=5.0u"eV",
                iter=0
            )
            s_walkers = [s_walker for _ in 1:3]
            liveset = LatticeGasWalkers(s_walkers, ham)

            e_type = typeof(s_walker.energy)
            
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=1000,              # mc_steps
                energy_perturbation=0.1,    # energy_perturbation
                fail_count=0,               # fail_count
                allowed_fail_count=100,     # allowed_fail_count
                random_seed=1234            # random_seed
            )
    
            @testset "MCRandomWalkMaxE" begin
                mc_routine = MCRandomWalkMaxE()
                iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                
                @test iter isa Union{Missing,Int}
                @test emax isa Union{Missing,e_type}
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test updated_params.fail_count >= 0
            end

            @testset "MCRandomWalkClone" begin
                mc_routine = MCRandomWalkClone()
                iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                
                @test iter isa Union{Missing,Int}
                @test emax isa Union{Missing,e_type}
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test updated_params.fail_count >= 0
            end
    
            @testset "MCNewSample" begin
                mc_routine = MCNewSample()
                iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                
                @test iter isa Union{Missing,Int}
                @test emax isa Union{Missing,e_type}
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test updated_params.fail_count >= 0
            end

            @testset "MCRejection" begin
                mc_routine = MCRejectionSampling()
                iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                
                @test iter isa Union{Missing,Int}
                @test emax isa Union{Missing,e_type}
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test updated_params.fail_count >= 0
            end
        end
    end


    @testset "adjust_step_size function tests" begin

        ns_params = NestedSamplingParameters(
                1000,    # mc_steps
                0.1,     # initial_step_size
                0.1,     # step_size
                0.01,    # step_size_lo
                1.0,     # step_size_up
                (0.5,0.5),  # acceptance ratio range
                0,       # fail_count
                100,     # allowed_fail_count
                1e-12,   # energy_perturbation
                1234     # random_seed
            )

        @testset "Step size increases" begin

            updated_params = deepcopy(ns_params)

            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.8)
            @test updated_params.step_size == ns_params.step_size * 1.1 # 10% increase
            
            # Test upper boundary
            updated_params.step_size = 0.8
            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.8)
            @test updated_params.step_size ≤ ns_params.step_size_up
        end
    
        @testset "Step size decreases" begin                
            
            updated_params = deepcopy(ns_params)

            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.2)
            @test updated_params.step_size == 0.1 * 0.9 # 10% decrease
            
            # Test lower boundary
            updated_params.step_size = 0.106
            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.2)
            @test updated_params.step_size ≥ ns_params.step_size_lo
        end
    
        @testset "No adjustment needed" begin
            
            updated_params = deepcopy(ns_params)

            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.5)
            @test updated_params.step_size == 0.1
            
            # Custom range
            updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.4, range=(0.3, 0.6))
            @test updated_params.step_size == 0.1
        end
    end


    @testset "nested_sampling_loop functions tests" begin
        @testset "AtomWalkers cases tests" begin
            # Setup test system
            box = [[10.0u"Å", 0u"Å", 0u"Å"],
                    [0u"Å", 10.0u"Å", 0u"Å"],
                    [0u"Å", 0u"Å", 10.0u"Å"]]
            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
            
            coor_list = [:H => [0.2, 0.5, 0.5],
                            :H => [0.4, 0.5, 0.5]]
            at = FastSystem(periodic_system(coor_list, box, fractional=true))
            
            walkers = [AtomWalker(at) for _ in 1:3]
            liveset = LJAtomWalkers(walkers, lj)
            
            ns_params = NestedSamplingParameters(
                1000,    # mc_steps
                0.1,     # initial_step_size
                0.1,     # step_size
                0.01,    # step_size_lo
                1.0,     # step_size_up
                (0.5,0.5),  # acceptance ratio range
                0,       # fail_count
                100,     # allowed_fail_count
                1e-12,   # energy_perturbation
                1234     # random_seed
            )
            
            save_strategy = SaveEveryN(
                df_filename = "test_df.csv",
                wk_filename = "test.traj.extxyz",
                ls_filename = "test.ls.extxyz",
                n_traj = 2,
                n_snap = 2,
                n_info = 2
            )

            @testset "Basic functionality" begin
                n_steps = 5
                ns_params_copy = deepcopy(ns_params)
                df, updated_liveset, updated_params = nested_sampling(
                    liveset, ns_params_copy, n_steps, MCRandomWalkMaxE(), save_strategy)
                
                @test df isa DataFrame
                @test names(df) == ["iter", "emax"]
                @test nrow(df) ≤ n_steps
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test eltype(df.iter) == Int
                @test eltype(df.emax) == Float64
            end

            @testset "Parallel Routine" begin
                n_steps = 5
                ns_params_copy = deepcopy(ns_params)
                df, updated_liveset, updated_params = nested_sampling(
                    liveset, ns_params_copy, n_steps, MCRandomWalkCloneParallel(), save_strategy)
                
                @test df isa DataFrame
                @test names(df) == ["iter", "emax"]
                @test nrow(df) ≤ n_steps
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test eltype(df.iter) == Int
                @test eltype(df.emax) == Float64
            end

            @testset "Failure handling" begin
                ns_params_copy = deepcopy(ns_params)
                ns_params_copy.allowed_fail_count = 1
                df, updated_liveset, updated_params = nested_sampling(
                    liveset, ns_params_copy, 10, MCRandomWalkMaxE(), save_strategy)
                
                @test updated_params.fail_count == 0
            end
    
            @testset "Data saving" begin
                ns_params_copy = deepcopy(ns_params)
                df, _, _ = nested_sampling(liveset, ns_params_copy, 4, MCRandomWalkMaxE(), save_strategy)
                
                @test isfile("test_df.csv")
                @test isfile("test.traj.extxyz")
                @test isfile("test.ls.extxyz")
                
                rm("test_df.csv", force=true)
                rm("test.traj.extxyz", force=true)
                rm("test.ls.extxyz", force=true)
            end

            @testset "Data saving with .arrow and SaveFreePartEveryN" begin
                ns_params_copy = deepcopy(ns_params)
                save_strategy = SaveFreePartEveryN(
                    df_filename = "test_df.arrow",
                    wk_filename = "test.traj.extxyz",
                    ls_filename = "test.ls.extxyz",
                    n_traj = 2,
                    n_snap = 2,
                    n_info = 2
                )
                df, _, _ = nested_sampling(liveset, ns_params_copy, 4, MCRandomWalkMaxE(), save_strategy)
                
                @test isfile("test_df.arrow")
                @test isfile("test.traj.extxyz")
                @test isfile("test.ls.extxyz")
                
                rm("test_df.arrow", force=true)
                rm("test.traj.extxyz", force=true)
                rm("test.ls.extxyz", force=true)
            end
        end

        @testset "LatticeGasWalkers cases tests" begin
            # Setup test system
            square_lattice = MLattice{1,SquareLattice}(
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(4, 4, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.5],
                components=:equal,
                adsorptions=:full
            )
            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            s_walker = LatticeWalker(square_lattice, energy=5.0u"eV", iter=0)
            s_walkers = [deepcopy(s_walker) for _ in 1:3]
            liveset = LatticeGasWalkers(s_walkers, ham)
            
            ns_params = LatticeNestedSamplingParameters(
                mc_steps=1000,
                energy_perturbation=0.1,
                fail_count=0,
                allowed_fail_count=100,
                random_seed=1234
            )
            save_strategy = SaveEveryN("test_df.csv", "test.traj", "test.ls", 2, 2, 2)

            @testset "Basic functionality" begin
                df, updated_liveset, updated_params = nested_sampling(
                    liveset, deepcopy(ns_params), 5, MCRandomWalkMaxE(), save_strategy)
                
                @test df isa DataFrame
                @test names(df) == ["iter", "emax"]
                @test eltype(df.emax) == Float64  # emax is stored as Float64
                @test length(updated_liveset.walkers) == length(liveset.walkers)
                @test all(walker -> walker isa LatticeWalker, updated_liveset.walkers)
            end

            @testset "Failure handling" begin
                fail_params = deepcopy(ns_params)
                fail_params.allowed_fail_count = 1
                df, _, updated_params = nested_sampling(
                    liveset, fail_params, 10, MCRandomWalkMaxE(), save_strategy)
                
                @test updated_params.fail_count == 0
                @test nrow(df) <= 10
            end

            @testset "Data saving" begin
                nested_sampling(liveset, deepcopy(ns_params), 4, MCRandomWalkMaxE(), save_strategy)
                
                @test isfile("test_df.csv")
                @test isfile("test.ls")
                
                rm("test_df.csv", force=true)
                rm("test.ls", force=true)
            end

            @testset "Walker properties" begin
                _, updated_liveset, _ = nested_sampling(
                    liveset, deepcopy(ns_params), 3, MCRandomWalkMaxE(), save_strategy)
                
                walker = updated_liveset.walkers[1]
                @test walker isa LatticeWalker
                @test typeof(walker.energy) <: Quantity
                @test unit(walker.energy) == u"eV"
                @test walker.configuration isa SLattice{SquareLattice}
            end

            # clean up
            if isfile("test_df.csv")
                rm("test_df.csv", force=true)
            end
            if isfile("test.traj")
                rm("test.traj", force=true)
            end

            if isfile("test.ls")
                rm("test.ls", force=true)
            end
        end
    end
end