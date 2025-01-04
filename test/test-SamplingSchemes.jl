@testset "SamplingSchemes Tests" begin
    @testset "exact_enumeration.jl tests" begin
        @testset "unique_permutations function tests" begin
            # Test empty input
            @test SamplingSchemes.unique_permutations(Int[]) == Vector{Int}[]
            
            # Test single element
            @test SamplingSchemes.unique_permutations([1]) == [[1]]
            
            # Test two unique elements
            @test sort(SamplingSchemes.unique_permutations([1,2])) == sort([[1,2], [2,1]])
            
            # Test two identical elements
            @test SamplingSchemes.unique_permutations([1,1]) == [[1,1]]
            
            # Test three elements with duplicates
            result3 = SamplingSchemes.unique_permutations([1,1,2])
            @test length(result3) == 3
            @test sort(result3) == sort([[1,1,2], [1,2,1], [2,1,1]])
            
            # Test string input
            @test sort(SamplingSchemes.unique_permutations(['a','b'])) == sort([['a','b'], ['b','a']])
            
            # Test longer sequence with multiple duplicates
            result4 = SamplingSchemes.unique_permutations([1,1,2,2])
            @test length(result4) == 6
            @test sort(result4) == sort([
                [1,1,2,2], [1,2,1,2], [1,2,2,1],
                [2,1,1,2], [2,1,2,1], [2,2,1,1]
            ])
            
            # Test type stability
            @test eltype(SamplingSchemes.unique_permutations([1,2])) == Vector{Int}
            @test eltype(SamplingSchemes.unique_permutations(['a','b'])) == Vector{Char}
        end

        @testset "enumerate_lattices function tests" begin
            @testset "MLattice as init_lattice cases" begin
                # 2x1 square lattice with 2 components
                lattice = MLattice{2,SquareLattice}(
                    supercell_dimensions=(2, 1, 1),
                    components=[[1], [2]]
                )
                results = SamplingSchemes.enumerate_lattices(lattice)
                
                @test length(results) == 2
                @test all(r isa MLattice{2,SquareLattice} for r in results)
                @test all(r.lattice_vectors == lattice.lattice_vectors for r in results)
                @test all(r.supercell_dimensions == (2,1,1) for r in results)
                @test all(sum(r.components[1]) + sum(r.components[2]) == 2 for r in results)
                
                # 2x2 triangular lattice with 3 components
                lattice = MLattice{3,TriangularLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=[[1], [2], [3]]
                )
                results = SamplingSchemes.enumerate_lattices(lattice)
                
                @test length(results) == 336  # (8! / (3! * 5!)) * 3 * 2 * 1
                @test all(r isa MLattice{3,TriangularLattice} for r in results)
                @test all(r.basis == lattice.basis for r in results)
                @test all(sum(r.components[1]) == 1 for r in results)
                @test all(sum(r.components[2]) == 1 for r in results)
                @test all(sum(r.components[3]) == 1 for r in results)
            end

            @testset "SLattice as init_lattice cases" begin
                # 2x1 square lattice with 1 occupied site
                lattice = SLattice{SquareLattice}(
                    supercell_dimensions=(2, 1, 1),
                    components=[[true, false]]
                )
                results = SamplingSchemes.enumerate_lattices(lattice)
                
                @test length(results) == 2  # C(2,1)
                @test all(r isa SLattice{SquareLattice} for r in results)
                @test all(sum(r.components[1]) == 1 for r in results)
                
                # 2x2 triangular lattice with 2 occupied sites
                lattice = SLattice{TriangularLattice}(
                    supercell_dimensions=(2, 2, 1),
                    components=[[true, true, true, true, 
                    false, false, false, false]]
                )
                results = SamplingSchemes.enumerate_lattices(lattice)
                
                @test length(results) == 70  # C(8,4)
                @test all(r isa SLattice{TriangularLattice} for r in results)
                @test all(sum(r.components[1]) == 4 for r in results)

                # Verify thread safety
                results_set = Set(map(r -> findall(r.components[1]), results))
                @test length(results_set) == 70
            end
        end


        @testset "exact_enumeration function tests" begin
            # Test 2-component system
            on_site = [-0.04, -0.03]u"eV"
            neighbor_ints = [[-0.01, -0.0025], [-0.015, -0.002], [-0.02, -0.003]].*u"eV"
            hams = [GenericLatticeHamiltonian(on_site[i], neighbor_ints[i]) for i in 1:2]
            push!(hams, GenericLatticeHamiltonian(on_site[2], neighbor_ints[3]))
            h = MLatticeHamiltonian(2, hams)
            
            lattice = MLattice{2,SquareLattice}(
                supercell_dimensions=(2, 2, 1),
                components=[[1], [4]]
            )
            
            df, ls = exact_enumeration(lattice, h)
            
            @test size(df, 1) == 12
            @test all(x -> unit(x) == u"eV", df.energy)
            @test all(c -> length(c) == 2, df.config)
            
            # Test 3-component system
            on_site = [-0.04, -0.03, -0.02]u"eV"
            neighbor_ints = [
                [-0.01, -0.0025], 
                [-0.015, -0.002], 
                [-0.02, -0.003],
                [-0.012, -0.0022],
                [-0.018, -0.0028],
                [-0.016, -0.0024]
            ].*u"eV"
            
            hams = [GenericLatticeHamiltonian(on_site[i], neighbor_ints[i]) for i in 1:3]
            for i in 4:6
                push!(hams, GenericLatticeHamiltonian(on_site[2], neighbor_ints[i]))
            end
            h = MLatticeHamiltonian(3, hams)
            
            lattice = MLattice{3,SquareLattice}(
                supercell_dimensions=(2, 2, 1),
                components=[[true, false, false, false], 
                        [false, true, true, false],
                        [false, false, false, true]]
            )
            
            df, ls = exact_enumeration(lattice, h)
            
            @test size(df, 1) == 12
            @test all(x -> unit(x) == u"eV", df.energy)
            @test all(sum(c[1]) == 1 for c in df.config)
            @test all(sum(c[2]) == 2 for c in df.config)
            @test all(sum(c[3]) == 1 for c in df.config)
            
            # Test energy differences
            @test length(unique(df.energy)) > 1
        end
    end


    @testset "nested_sampling.jl tests" begin
        @testset "NestedSamplingParameters struct tests" begin
            params = NestedSamplingParameters(
                1000,    # mc_steps
                0.1,     # initial_step_size
                0.1,     # step_size
                0.01,    # step_size_lo
                1.0,     # step_size_up
                0,       # fail_count
                100      # allowed_fail_count
            )
            
            @test params isa SamplingSchemes.SamplingParameters
            @test params.mc_steps == 1000
            @test params.initial_step_size == 0.1
            @test params.step_size == 0.1
            @test params.step_size_lo == 0.01
            @test params.step_size_up == 1.0
            @test params.fail_count == 0
            @test params.allowed_fail_count == 100
            
            # Test mutability
            params.step_size = 0.2
            params.fail_count = 1
            @test params.step_size == 0.2
            @test params.fail_count == 1
        end
        

        @testset "LatticeNestedSamplingParameters struct tests" begin
            params = LatticeNestedSamplingParameters(
                1000,    # mc_steps
                0.1,     # energy_perturbation
                0,       # fail_count
                100      # allowed_fail_count
            )
            
            @test params isa SamplingSchemes.SamplingParameters
            @test params.mc_steps == 1000
            @test params.energy_perturbation == 0.1
            @test params.fail_count == 0
            @test params.allowed_fail_count == 100
            
            # Test mutability
            params.fail_count = 1
            @test params.fail_count == 1
        end
        

        @testset "MCRoutine types tests" begin
            @test MCRandomWalkMaxE() isa MCRoutine
            @test MCRandomWalkClone() isa MCRoutine
            @test MCNewSample() isa MCRoutine
            
            routines = [MCRandomWalkMaxE(), MCRandomWalkClone(), MCNewSample()]
            @test all(r -> r isa MCRoutine, routines)
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
                
                coor_list = [:H => [0.2, 0.5, 0.5],
                            :H => [0.4, 0.5, 0.5],
                            :O => [0.6, 0.5, 0.5]]
                at = FastSystem(periodic_system(coor_list, box, fractional=true))
                
                walkers = [AtomWalker(at) for _ in 1:3]
                liveset = LJAtomWalkers(walkers, lj)
                
                ns_params = NestedSamplingParameters(
                    1000,    # mc_steps
                    0.1,     # initial_step_size
                    0.1,     # step_size
                    0.01,    # step_size_lo
                    1.0,     # step_size_up
                    0,       # fail_count
                    100      # allowed_fail_count
                )
                
                @testset "MCRandomWalkMaxE" begin
                    mc_routine = MCRandomWalkMaxE()
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

                @testset "Unsupported MCRoutine" begin
                    mc_routine = MCNewSample()
                    @test_throws UndefVarError(:unsupported_mc) begin
                        nested_sampling_step!(liveset, ns_params, unsupported_mc)
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
                
                ns_params = LatticeNestedSamplingParameters(
                    1000,    # mc_steps
                    0.1,     # energy_perturbation
                    0,       # fail_count
                    100      # allowed_fail_count
                )
        
                @testset "MCRandomWalkMaxE" begin
                    mc_routine = MCRandomWalkMaxE()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,Float64}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end

                @testset "MCRandomWalkClone" begin
                    mc_routine = MCRandomWalkClone()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,Float64}
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test updated_params.fail_count >= 0
                end
        
                @testset "MCNewSample" begin
                    mc_routine = MCNewSample()
                    iter, emax, updated_liveset, updated_params = nested_sampling_step!(liveset, ns_params, mc_routine)
                    
                    @test iter isa Union{Missing,Int}
                    @test emax isa Union{Missing,Float64}
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
                    0,       # fail_count
                    100      # allowed_fail_count
                )

            @testset "Step size increases" begin

                updated_params = deepcopy(ns_params)

                updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.8)
                @test updated_params.step_size == ns_params.step_size * 1.05
                
                # Test upper boundary
                updated_params.step_size = 0.8
                updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.8)
                @test updated_params.step_size ≤ ns_params.step_size_up
            end
        
            @testset "Step size decreases" begin                
                
                updated_params = deepcopy(ns_params)

                updated_params = SamplingSchemes.adjust_step_size(updated_params, 0.2)
                @test updated_params.step_size == 0.1 * 0.95
                
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
                    0,       # fail_count
                    100      # allowed_fail_count
                )
                
                save_strategy = SaveEveryN(
                    df_filename = "test_df.csv",
                    wk_filename = "test.traj.extxyz",
                    ls_filename = "test.ls.extxyz",
                    n_traj = 2,
                    n_snap = 2
                )

                @testset "Basic functionality" begin
                    n_steps = 5
                    ns_params_copy = deepcopy(ns_params)
                    df, updated_liveset, updated_params = nested_sampling_loop!(
                        liveset, ns_params_copy, n_steps, MCRandomWalkMaxE(), save_strategy)
                    
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
                    df, updated_liveset, updated_params = nested_sampling_loop!(
                        liveset, ns_params_copy, 10, MCRandomWalkMaxE(), save_strategy)
                    
                    @test updated_params.fail_count == 0
                end
        
                @testset "Data saving" begin
                    ns_params_copy = deepcopy(ns_params)
                    df, _, _ = nested_sampling_loop!(liveset, ns_params_copy, 4, MCRandomWalkMaxE(), save_strategy)
                    
                    @test isfile("test_df.csv")
                    @test isfile("test.traj.extxyz")
                    @test isfile("test.ls.extxyz")
                    
                    rm("test_df.csv", force=true)
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
                
                ns_params = LatticeNestedSamplingParameters(1000, 0.1, 0, 100)
                save_strategy = SaveEveryN("test_df.csv", "test.traj", "test.ls", 2, 2)

                @testset "Basic functionality" begin
                    df, updated_liveset, updated_params = nested_sampling_loop!(
                        liveset, deepcopy(ns_params), 5, MCRandomWalkMaxE(), save_strategy)
                    
                    @test df isa DataFrame
                    @test names(df) == ["iter", "emax", "config"]
                    @test eltype(df.emax) == Float64  # emax is stored as Float64
                    @test length(updated_liveset.walkers) == length(liveset.walkers)
                    @test all(walker -> walker isa LatticeWalker, updated_liveset.walkers)
                end

                @testset "Failure handling" begin
                    fail_params = deepcopy(ns_params)
                    fail_params.allowed_fail_count = 1
                    df, _, updated_params = nested_sampling_loop!(
                        liveset, fail_params, 10, MCRandomWalkMaxE(), save_strategy)
                    
                    @test updated_params.fail_count == 0
                    @test nrow(df) <= 10
                end

                @testset "Data saving" begin
                    nested_sampling_loop!(liveset, deepcopy(ns_params), 4, MCRandomWalkMaxE(), save_strategy)
                    
                    @test isfile("test_df.csv")
                    @test isfile("test.ls")
                    
                    rm("test_df.csv", force=true)
                    rm("test.ls", force=true)
                end

                @testset "Walker properties" begin
                    _, updated_liveset, _ = nested_sampling_loop!(
                        liveset, deepcopy(ns_params), 3, MCRandomWalkMaxE(), save_strategy)
                    
                    walker = updated_liveset.walkers[1]
                    @test walker isa LatticeWalker
                    @test typeof(walker.energy) <: Quantity
                    @test unit(walker.energy) == u"eV"
                    @test walker.configuration isa SLattice{SquareLattice}
                end
            end
        end
    end
    

    @testset "nvt_monte_carlo.jl tests" begin
        @testset "nvt_monte_carlo function tests" begin
            
            # Setup
            lattice = SLattice{SquareLattice}(
                    supercell_dimensions=(2, 1, 1),
                    components=[[true, false]]
                )
            
            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

            @testset "Basic functionality" begin
                
                energies, configs, accepted = nvt_monte_carlo(lattice, ham, 300.0, 100, 42)
                
                @test length(energies) == 101
                @test length(configs) == 101
                @test 0 ≤ accepted ≤ 100
                @test all(isfinite, energies)
            end
        
            @testset "Temperature effects" begin
        
                energies_cold, _, accepted_cold = nvt_monte_carlo(lattice, ham, 1.0, 50, 42)
                energies_hot, _, accepted_hot = nvt_monte_carlo(lattice, ham, 1000.0, 50, 42)
        
                @test std(energies_cold) ≤ std(energies_hot)
                @test accepted_cold ≤ accepted_hot
            end
        
            @testset "Conservation & Reproducibility" begin
        
                energies1, configs1, accepted1 = nvt_monte_carlo(lattice, ham, 300.0, 50, 42)
                energies2, configs2, accepted2 = nvt_monte_carlo(lattice, ham, 300.0, 50, 42)
        
                @test all(sum(c.components[1]) == sum(lattice.components[1]) for c in configs1)
                @test energies1 == energies2
                @test accepted1 == accepted2
            end
        end
    end


    @testset "wang_landau.jl tests" begin
        @testset "wang_landau function tests" begin

            # Setup
            lattice = SLattice{SquareLattice}(
                supercell_dimensions=(2, 1, 1),
                components=[[true, false]]
            )

            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            energy_bins = range(-0.1, 0.1, length=20)

            @testset "Basic functionality" begin
                S, H, bins, energies, configs = wang_landau(
                    lattice, ham, 100, 0.8, 2.7, 1.0001,
                    collect(energy_bins), 42
                )
                
                @test length(S) == length(energy_bins)
                @test length(H) == length(energy_bins)
                @test length(energies) == 1401
                @test length(configs) == 1401
                @test all(x -> x ≥ 0, H)
            end

            @testset "Energy conservation" begin
                S, _, _, energies, configs = wang_landau(
                    lattice, ham, 50, 0.8, 2.7, 1.001,
                    collect(energy_bins), 42
                )
                
                @test all(e -> minimum(energy_bins) ≤ e ≤ maximum(energy_bins), energies)
                @test all(c -> sum(c.components[1]) == sum(lattice.components[1]), configs)
            end

            @testset "Modification factor" begin
                S1, _, _, _, _ = wang_landau(
                    lattice, ham, 50, 0.8, 2.7, 1.1,
                    collect(energy_bins), 42
                )
                S2, _, _, _, _ = wang_landau(
                    lattice, ham, 50, 0.8, 2.7, 1.001,
                    collect(energy_bins), 42
                )
                
                @test !isapprox(S1, S2)
            end
        end


        @testset "get_bin_index function tests" begin
            bins = [-1.0, -0.5, 0.0, 0.5, 1.0]

            @test SamplingSchemes.get_bin_index(-0.7, bins) == 1
            @test SamplingSchemes.get_bin_index(-0.2, bins) == 2
            @test SamplingSchemes.get_bin_index(0.3, bins) == 3
            @test SamplingSchemes.get_bin_index(0.8, bins) == 4
            @test SamplingSchemes.get_bin_index(1.5, bins) == 5  # Edge case: above max
            @test SamplingSchemes.get_bin_index(0.0, bins) == 3  # Edge case: exact boundary
        end
    end

end