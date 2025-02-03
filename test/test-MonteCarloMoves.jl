@testset "Monte Carlo Moves tests" begin
    @testset "random_walks.jl tests" begin
        @testset "single_atom_random_walk function tests" begin

            # Test basic functionality
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 2.0, 3.0] .* u"Å")
            step_size = 0.5
            new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
            
            @test new_pos ≠ pos
            
            # Test step size constraints
            for _ in 1:10
                pos = SVector{3, typeof(0.0u"Å")}([0.0, 0.0, 0.0] .* u"Å")
                new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
                
                # Check that no component moved more than step_size
                @test all(abs.(new_pos - pos) .<= step_size * u"Å")
            end
            
            # Test type stability
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 1.0, 1.0] .* u"Å")
            @test typeof(MonteCarloMoves.single_atom_random_walk!(pos, step_size)) == typeof(pos)
        end


        @testset "MC_random_walk functions tests" begin
            @testset "AtomWalker MC_random_walk tests" begin
                
                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0u"eV"
                
                # Create a simple test system
                at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])    # Example values
                atw = AtomWalker(at;freeze_species=[:H],merge_same_species=true)
                lj = LJParameters()  # Create test LJ parameters
                
                # Run MC simulation
                accept_this_walker, accept_rate, updated_atw = MC_random_walk!(
                    n_steps, atw, lj, step_size, emax
                )
                
                # Tests
                @test typeof(accept_this_walker) == Bool
                @test 0.0 <= accept_rate <= 1.0
                @test typeof(updated_atw) == typeof(atw)
                @test updated_atw.energy <= emax
            end

            @testset "LatticeWalker MC_random_walk tests" begin

                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0             # Float64
                
                # Create a simple square lattice system for test
                square_lattice = MLattice{1,SquareLattice}(
                    lattice_constant=1.0,                   # Unit cell size
                    basis=[(0.0, 0.0, 0.0)],                # Single atom per unit cell
                    supercell_dimensions=(4, 4, 1),         # 4x4x1 supercell
                    periodicity=(true, true, false),        # Periodic in x,y but not z
                    cutoff_radii=[1.1, 1.5],                # Neighbor cutoff distances
                    components=:equal,                      # Split into equal components
                    adsorptions=:full                       # All sites are adsorption sites
                )

                s_walker = LatticeWalker(
                    square_lattice,
                    energy=5.0u"eV",
                    iter=0
                )

                # Create a simple triangular lattice system for test
                triangular_lattice = MLattice{1,TriangularLattice}(
                    lattice_constant=1.0,
                    basis=[(0.0, 0.0, 0.0), (1/2, sqrt(3)/2, 0.0)],  # Two atoms per unit cell
                    supercell_dimensions=(4, 2, 1),
                    periodicity=(true, true, false),
                    cutoff_radii=[1.1, 1.5],
                    components=:equal,
                    adsorptions=:full
                )

                t_walker = LatticeWalker(
                    triangular_lattice,
                    energy=5.0u"eV",
                    iter=0
                )

                # Create a simple Hamiltonian for testing
                ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

                # Test basic functionality
                @testset "Basic MC walk tests" begin
                    # Store initial state
                    initial_config = deepcopy(s_walker.configuration)
                    initial_energy = s_walker.energy
                    
                    # Perform MC walk
                    accepted, rate, updated_walker = MC_random_walk!(n_steps, s_walker, ham, emax)            # !!!Error: type MLattice has no field occupations
                    
                    # Basic checks
                    @test typeof(accepted) == Bool
                    @test 0.0 <= rate <= 1.0
                    @test updated_walker isa LatticeWalker
                end
            end
        end

        
        @testset "MC_new_sample function tests" begin

            # Create base test lattices
            square_lattice = MLattice{1,SquareLattice}(
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(4, 4, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.5],
                components=:equal,
                adsorptions=:full
            )

            lattice = LatticeWalker(
                square_lattice,
                energy=5.0u"eV",
                iter=0
            )
            
            h = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")

            @testset "Basic functionality and energy constraints tests" begin

                # Test with normal emax
                emax = 2.0
                result, updated_lattice = MC_new_sample!(lattice, h, emax)
                @test isa(result, Bool)
                @test isa(updated_lattice, LatticeWalker)
                @test unit(updated_lattice.energy) == u"eV"
                
                # Test with very low emax to force rejection
                result_low, _ = MC_new_sample!(lattice, h, -10.0)
                @test result_low == false
            end

            @testset "Energy perturbation and configuration updates tests" begin

                emax = 2.0
                # Store initial configuration
                original_config = deepcopy(lattice.configuration)
                
                # Test with perturbation
                result, updated_lattice = MC_new_sample!(lattice, h, emax, energy_perturb=0.1)
                @test isa(result, Bool)
                
                if result
                    @test updated_lattice.configuration ≠ original_config
                else
                    @test updated_lattice.configuration == original_config
                end
            end

            @testset "Different lattice type tests" begin

                # Test with triangular lattice
                tri_lattice = MLattice{1,TriangularLattice}(
                    lattice_constant=1.0,
                    basis=[(0.0, 0.0, 0.0), (1/2, sqrt(3)/2, 0.0)],
                    supercell_dimensions=(4, 2, 1),
                    periodicity=(true, true, false),
                    cutoff_radii=[1.1, 1.5],
                    components=:equal,
                    adsorptions=:full
                )
        
                t_walker = LatticeWalker(tri_lattice, energy=5.0u"eV", iter=0)
                result, _ = MC_new_sample!(t_walker, h, 2.0)
                @test isa(result, Bool)
            end

        end


        @testset "generate_random_new_lattice_sample function tests" begin
            @testset "Single-component lattice cases" begin

                # Define basic lattice parameters
                lattice_vectors = [1.0 0.0 0.0;
                                    0.0 1.0 0.0;
                                    0.0 0.0 1.0]
                basis = [(0.0, 0.0, 0.0)]
                supercell_dims = (4, 4, 1)
                periodicity = (true, true, false)
                cutoff_radii = [1.1]
                
                # Create components for single-component system (4 sites total)
                components = [fill(true, 4)]  # One component, all sites initially occupied
                adsorptions = fill(true, 4)   # All sites are adsorption sites
                
                square_lattice = MLattice{1,SquareLattice}(
                    lattice_vectors,
                    basis,
                    supercell_dims,
                    periodicity,
                    cutoff_radii,
                    components,
                    adsorptions
                )
                
                # Store initial occupancy
                initial_occupied = sum(square_lattice.components[1])
                
                # Generate new sample
                generate_random_new_lattice_sample!(square_lattice)
                
                # Test preservation of number of occupied sites
                @test sum(square_lattice.components[1]) == initial_occupied
                @test length(square_lattice.components) == 1  # Single component
            end
        
            @testset "Multi-component lattice cases" begin

                # Define same basic lattice parameters
                lattice_vectors = [1.0 0.0 0.0;
                                    0.0 1.0 0.0;
                                    0.0 0.0 1.0]
                basis = [(0.0, 0.0, 0.0)]
                supercell_dims = (2, 2, 1)
                periodicity = (true, true, false)
                cutoff_radii = [1.1, 1.5]
                
                # Create components for two-component system (4 sites total)
                components = [
                    [true, true, false, false],  # First component occupies first 2 sites
                    [false, false, true, true]   # Second component occupies last 2 sites
                ]
                adsorptions = fill(true, 4)
                
                two_comp_lattice = MLattice{2,SquareLattice}(
                    lattice_vectors,
                    basis,
                    supercell_dims,
                    periodicity,
                    cutoff_radii,
                    components,
                    adsorptions
                )
                
                # Store initial occupancies
                initial_occupancy = occupied_site_count(two_comp_lattice)
                
                # Generate new sample
                generate_random_new_lattice_sample!(two_comp_lattice)
                
                # Test occupancy preservation for each component
                final_occupancy = occupied_site_count(two_comp_lattice)
                @test initial_occupancy == final_occupancy
                
                # Test no site is occupied by multiple components
                total_sites = num_sites(two_comp_lattice)
                for site in 1:total_sites
                    @test sum(comp[site] for comp in two_comp_lattice.components) ≤ 1
                end
            end
        
            @testset "Empty occupied and fully occupied cases" begin

                lattice_vectors = [1.0 0.0 0.0;
                                  0.0 1.0 0.0;
                                  0.0 0.0 1.0]
                basis = [(0.0, 0.0, 0.0)]
                supercell_dims = (2, 2, 1)
                periodicity = (true, true, false)
                cutoff_radii = [1.1]
                
                # Empty lattice
                empty_components = [fill(false, 4)]  # One component, no sites occupied
                adsorptions = fill(true, 4)
                
                empty_lattice = MLattice{1,SquareLattice}(
                    lattice_vectors,
                    basis,
                    supercell_dims,
                    periodicity,
                    cutoff_radii,
                    empty_components,
                    adsorptions
                )
                
                generate_random_new_lattice_sample!(empty_lattice)
                @test all(.!empty_lattice.components[1])  # All sites should remain unoccupied
                
                # Fully occupied lattice
                full_components = [fill(true, 4)]  # One component, all sites occupied
                
                full_lattice = MLattice{1,SquareLattice}(
                    lattice_vectors,
                    basis,
                    supercell_dims,
                    periodicity,
                    cutoff_radii,
                    full_components,
                    adsorptions
                )
                
                generate_random_new_lattice_sample!(full_lattice)
                @test sum(full_lattice.components[1]) == 4  # Should maintain full occupancy
            end
        end
    end


    @testset "helpers.jl tests" begin

        # Set up basic periodic box
        box_size = 10.0u"Å"
        box = [[box_size, 0u"Å", 0u"Å"],
               [0u"Å", box_size, 0u"Å"],
               [0u"Å", 0u"Å", box_size]]

        @testset "periodic_boundary_wrap function tests" begin

            # Single-atom system
            coor_list = [:H => [0.5, 0.5, 0.5] .* box_size]
            
            # Test periodic wrapping
            pbc = (true, true, true)
            periodic_sys = atomic_system(coor_list, box, pbc)
            pos = SVector{3,typeof(1.0u"Å")}([12.0, 15.0, -2.0] .* u"Å")
            wrapped_pos = periodic_boundary_wrap!(pos, periodic_sys)
            @test wrapped_pos ≈ SVector{3,typeof(1.0u"Å")}([2.0, 5.0, 8.0] .* u"Å")

            # Test Dirichlet reflection
            pbc = (false, false, false)
            dirichlet_sys = atomic_system(coor_list, box, pbc)
            pos = SVector{3,typeof(1.0u"Å")}([12.0, -2.0, 15.0] .* u"Å")
            wrapped_pos = periodic_boundary_wrap!(pos, dirichlet_sys)
            @test wrapped_pos ≈ SVector{3,typeof(1.0u"Å")}([8.0, 2.0, 5.0] .* u"Å")

            # Test mixed boundaries
            pbc = (true, false, true)
            mixed_sys = atomic_system(coor_list, box, pbc)
            pos = SVector{3,typeof(1.0u"Å")}([12.0, 12.0, -2.0] .* u"Å")
            wrapped_pos = periodic_boundary_wrap!(pos, mixed_sys)
            @test wrapped_pos ≈ SVector{3,typeof(1.0u"Å")}([2.0, 8.0, 8.0] .* u"Å")
        end

            
        @testset "free_par_index function tests" begin

            # Multi-atom system with different types
            coor_list = [:H => [0.2, 0.3, 0.5], 
                            :H => [0.8, 0.3, 0.5],
                            :O => [0.5, 0.2, 0.3], 
                            :O => [0.7, 0.8, 0.3]]
            at = FastSystem(periodic_system(coor_list, box, fractional=true))
            
            # Test with different freezing configurations
            @test length(MonteCarloMoves.free_par_index(AtomWalker(at, freeze_species=[:O]))) == 2
            @test isempty(MonteCarloMoves.free_par_index(AtomWalker(at, freeze_species=[:H,:O])))
            @test length(MonteCarloMoves.free_par_index(AtomWalker(at, freeze_species=Symbol[]))) == 4
        end
    

        @testset "mean_sq_displacement function tests" begin

            # Test periodic boundary crossing
            sys1 = FastSystem(periodic_system([:H => [0.9, 0.9, 0.9]], box, fractional=true))
            sys2 = FastSystem(periodic_system([:H => [0.1, 0.1, 0.1]], box, fractional=true))
            @test MonteCarloMoves.mean_sq_displacement(AtomWalker(sys2), AtomWalker(sys1)) ≈ 12.0u"Å^2"
    
            # Test normal displacement
            sys3 = FastSystem(periodic_system([:H => [0.2, 0.2, 0.2]], box, fractional=true))
            sys4 = FastSystem(periodic_system([:H => [0.3, 0.3, 0.3]], box, fractional=true))
            @test MonteCarloMoves.mean_sq_displacement(AtomWalker(sys4), AtomWalker(sys3)) ≈ 3.0u"Å^2"
    
            # Test frozen particle error
            sys_frozen = FastSystem(periodic_system([:H => [0.5, 0.5, 0.5]], box, fractional=true))
            atw_frozen = AtomWalker(sys_frozen, freeze_species=[:H])
            atw_moved = deepcopy(atw_frozen)
            atw_moved.configuration.position[1] += [1.0, 1.0, 1.0]u"Å"
            @test_throws UndefVarError mean_sq_displacement(atw_moved, atw_frozen)
        end
    end

end

