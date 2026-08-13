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

        @testset "single_atom_random_walk 2D/1D tests" begin

            # Test basic functionality
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 2.0, 3.0] .* u"Å")
            step_size = 0.5
            new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size, [1,2])
            
            @test new_pos ≠ pos
            @test new_pos[3] == pos[3] # Check that z-coordinate is unchanged
            
            # Test step size constraints
            for _ in 1:10
                pos = SVector{3, typeof(0.0u"Å")}([0.0, 0.0, 0.0] .* u"Å")
                new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size,[1,2])
                
                # Check that no component moved more than step_size
                @test all(abs.(new_pos - pos) .<= step_size * u"Å")
                @test new_pos[3] == pos[3] # Check that z-coordinate is unchanged
            end
        
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

            @testset "AtomWalker MC_random_walk_2D! tests" begin
                
                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0u"eV"
                
                # Create a simple test system
                at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[2,1,3,4,5,6];particle_types=[:H,:O,:H,:Fe,:Au,:Cl])    # Example values
                atw = AtomWalker(at;freeze_species=[:H],merge_same_species=true)
                lj = LJParameters()  # Create test LJ parameters
                
                # Run MC simulation
                accept_this_walker, accept_rate, updated_atw = MC_random_walk_2D!(
                    n_steps, atw, lj, step_size, emax
                )
                
                # Tests
                @test typeof(accept_this_walker) == Bool
                @test 0.0 <= accept_rate <= 1.0
                @test typeof(updated_atw) == typeof(atw)
                @test updated_atw.energy <= emax

                @test updated_atw.configuration.position[3] == atw.configuration.position[3] # Check that z-coordinate is unchanged
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

        @testset "MC_rejection_sampling function tests" begin

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
                result, updated_lattice = MC_rejection_sampling!(lattice, h, emax)
                @test isa(result, Bool)
                @test isa(updated_lattice, LatticeWalker)
                @test unit(updated_lattice.energy) == u"eV"
                
                # Test with very low emax to force rejection
                result_low, _ = MC_rejection_sampling!(lattice, h, -10.0)
                @test result_low == false
            end

            @testset "Different lattice type tests" begin

                # Test with triangular lattice
                tri_lattice = MLattice{1,TriangularLattice}(
                    lattice_constant=1.0,
                    basis=[(0.0, 0.0, 0.0), (1/2, sqrt(3)/2, 0.0)],
                    supercell_dimensions=(4
                    , 2, 1),
                    periodicity=(true, true, false),
                    cutoff_radii=[1.1, 1.5],
                    components=:equal,
                    adsorptions=:full
                )

                t_walker = LatticeWalker(tri_lattice, energy=5.0u"eV", iter=0)
                result, _ = MC_rejection_sampling!(t_walker, h, 2.0)
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

        @testset "lattice_biased_sites tests" begin
            biased_lattice() = MLattice{1,SquareLattice}(
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(4, 4, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.5],   # shell 1: 4 NN at d=1; shell 2: 4 diagonals at d=sqrt(2)
                components=:equal,
                adsorptions=:full)

            @testset "empty and full lattices" begin
                lat = biased_lattice()
                lat.components[1] .= false
                @test sort(lattice_biased_sites(lat; predicate=:cavity)) == collect(1:16)
                @test isempty(lattice_biased_sites(lat; predicate=:contact))
                lat.components[1] .= true
                @test isempty(lattice_biased_sites(lat; predicate=:cavity))
                @test isempty(lattice_biased_sites(lat; predicate=:contact))
            end

            @testset "single particle shell counts" begin
                lat = biased_lattice()
                lat.components[1] .= false
                lat.components[1][6] = true   # torus: every site equivalent
                # shells=1: contacts are exactly the 4 NN; cavities 16-1-4 = 11
                contact1 = lattice_biased_sites(lat; predicate=:contact, shells=1)
                cavity1 = lattice_biased_sites(lat; predicate=:cavity, shells=1)
                @test sort(contact1) == sort(lat.neighbors[6][1])
                @test length(contact1) == 4
                @test length(cavity1) == 11
                # shells=2: 4 NN + 4 diagonals; cavities 16-1-8 = 7
                contact2 = lattice_biased_sites(lat; predicate=:contact, shells=2)
                cavity2 = lattice_biased_sites(lat; predicate=:cavity, shells=2)
                @test sort(contact2) == sort(vcat(lat.neighbors[6][1], lat.neighbors[6][2]))
                @test length(contact2) == 8
                @test length(cavity2) == 7
            end

            @testset "contact and cavity partition the empty sites" begin
                lat = biased_lattice()
                lat.components[1] .= false
                lat.components[1][[1, 4, 6, 7, 10, 13, 16]] .= true  # fixed irregular config
                for shells in (1, 2)
                    contact = lattice_biased_sites(lat; predicate=:contact, shells=shells)
                    cavity = lattice_biased_sites(lat; predicate=:cavity, shells=shells)
                    @test isempty(intersect(contact, cavity))
                    @test sort(vcat(contact, cavity)) == findall(.!lat.components[1])
                end
            end

            @testset "image-multiplicity neighbor lists" begin
                # Faithful 4x4 cell: multiplicity lists equal min-image lists,
                # so the counts must be unchanged
                lat_m = MLattice{1,SquareLattice}(
                    lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
                    supercell_dimensions=(4, 4, 1), periodicity=(true, true, false),
                    cutoff_radii=[1.1, 1.5], components=[[false for _ in 1:16]],
                    adsorptions=:full, image_multiplicity=true)
                lat_m.components[1][6] = true
                @test length(lattice_biased_sites(lat_m; predicate=:contact)) == 4
                @test length(lattice_biased_sites(lat_m; predicate=:cavity)) == 11
                # 2x2 torus: each site's 4 NN collapse onto 2 distinct sites
                # (double entries); the predicates read occupancy only, so
                # duplicated neighbor entries are harmless
                tiny = MLattice{1,SquareLattice}(
                    lattice_constant=1.0, basis=[(0.0, 0.0, 0.0)],
                    supercell_dimensions=(2, 2, 1), periodicity=(true, true, false),
                    cutoff_radii=[1.1], components=[[false for _ in 1:4]],
                    adsorptions=:full, image_multiplicity=true)
                tiny.components[1][1] = true
                @test sort(lattice_biased_sites(tiny; predicate=:contact)) == [2, 3]
                @test lattice_biased_sites(tiny; predicate=:cavity) == [4]
            end

            @testset "argument validation" begin
                lat = biased_lattice()
                @test_throws ArgumentError lattice_biased_sites(lat; predicate=:nonsense)
                @test_throws ArgumentError lattice_biased_sites(lat; shells=0)
                @test_throws ArgumentError lattice_biased_sites(lat; shells=3)  # lattice has 2 shells
            end
        end
    end

    @testset "Monte Carlo swap moves tests" begin
        @testset "two_atoms_swap! function tests" begin
            
            # Create a simple test system
            at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[1,1,1,1,1];particle_types=[:H,:O,:Fe,:Au,:Cl])    # Example values
            atw = AtomWalker(at)
            lj = LJParameters()  # Create test LJ parameters
            
            # Run MC simulation
            atw_new = MonteCarloMoves.two_atoms_swap!(deepcopy(atw), 1, 2)
            
            # Tests
            @test atw_new.configuration.position[1] == atw.configuration.position[2]
            @test atw_new.configuration.position[2] == atw.configuration.position[1]
        end

        @testset "MC_random_swap! function tests" begin
            
            # Create a simple test system
            at = FreeBirdIO.generate_multi_type_random_starting_config(10.0,[1,1,1,1,1];particle_types=[:H,:O,:Fe,:Au,:Cl])    # Example values
            atw = AtomWalker(at;freeze_species=[:Fe,:Au,:Cl],merge_same_species=true) # only H and O are free
            lj = LJParameters()  # Create test LJ parameters
            
            # Run MC simulation
            accepted, rate, atw_new = MonteCarloMoves.MC_random_swap!(1, deepcopy(atw), lj, Inf*u"eV") # Inf energy limit, so always accept
            
            # Tests
            @test atw_new.configuration.position[1] == atw.configuration.position[2]
            @test atw_new.configuration.position[2] == atw.configuration.position[1]
            @test atw_new.configuration.position[3] == atw.configuration.position[3] # Check if other atoms are unchanged

            @test typeof(accepted) == Bool
            @test rate == 1.0
            @test accepted == true

            accepted, rate, atw_new = MonteCarloMoves.MC_random_swap!(1, deepcopy(atw), lj, -Inf*u"eV") # -Inf energy limit, so always reject

            @test atw_new.configuration.position[1] == atw.configuration.position[1] # Check if positions are unchanged
            @test atw_new.configuration.position[2] == atw.configuration.position[2] # Check if positions are unchanged
            @test atw_new.configuration.position[3] == atw.configuration.position[3] # Check if other atoms are unchanged

            @test typeof(accepted) == Bool
            @test rate == 0.0
            @test accepted == false
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

        @testset "free_component_index function tests" begin

            # Multi-atom system with different types
            coor_list = [:H => [0.2, 0.3, 0.5], 
                            :H => [0.8, 0.3, 0.5],
                            :O => [0.5, 0.2, 0.3], 
                            :O => [0.7, 0.8, 0.3]]
            at = FastSystem(periodic_system(coor_list, box, fractional=true))
            
            # Test with different freezing configurations
            @test length(MonteCarloMoves.free_component_index(AtomWalker(at, freeze_species=[:O]))) == 1
            @test isempty(MonteCarloMoves.free_component_index(AtomWalker(at, freeze_species=[:H,:O])))
            @test length(MonteCarloMoves.free_component_index(AtomWalker(at, freeze_species=Symbol[]))) == 2
        end
    end

    # Calibration ledger for the atomistic grand-canonical kernel testsets below
    # (protocol: burn 2e4 steps, then 2e5 steps sampling N every 10; per-attempt
    # acceptance references are attempt-conditioned, i.e. deletion guard skips at
    # N = 0 excluded, matching the kernel's counting contract).
    # Corner z0V = 0.5, seeds {1,2,3,42421}: max devs mean 3.85e-3, var 4.72e-3,
    #   P(0) 3.97e-3, A_ins 3.24e-3; A_del is structurally exact 1.0 at this corner.
    # Corner z0V = 8.0, seeds {1,2,3,42422}: max devs mean 3.31e-2, var 3.27e-1,
    #   P(0) 1.35e-4, A_ins 1.59e-3, A_del(conditioned) 8.7e-3 (gate 0.027 >= 3x).
    # Corner z0V = 32.0, seeds {1,2,3,42423}: max devs mean 3.41e-1, var 1.11,
    #   A_ins 3.62e-3, A_del(conditioned) 5.9e-3.
    # Two-state n_max = 1 at z0V = 0.6, seeds {1,2,3,42424}: max dev f1 4.26e-3.
    # The zero-interaction digit pin (seed 424242) is architecture-exact (every
    # energy is exactly 0.0 eV, so no rounding enters the trajectory); the
    # interacting pin (seed 424243) gates its energy at rtol 1e-12 to tolerate
    # compiler-level rounding drift along the same trajectory.
    # Gates ship at >= 3x the max deviation per stat per corner. Off-by-one foil
    # sensitivity (exact birth-death mean at z0V = 0.5): shift +0.1778, so the
    # mean gate 0.012 detects it with ~15x headroom (asserted in-test).
    @testset "atomistic grand-canonical kernel tests" begin
        using Random
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
        pbc = (true, true, true)
        seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
        mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                               empty(position(seed_at, :)), empty(species(seed_at, :)),
                               empty(mass(seed_at, :)))
        mkwalker() = (w = AtomWalker{1}(mkempty()); w.energy = 0.0u"eV"; w)
        lj0 = LJParameters(epsilon=0.0)
        emax_inf = Inf * u"eV"

        logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
        poisson_pmf(lam, n) = exp(n * log(lam) - lam - logfact(n))
        function poisson_refs(lam; nmax=max(80, ceil(Int, lam + 12 * sqrt(lam))))
            p = [poisson_pmf(lam, n) for n in 0:nmax]
            p ./= sum(p)
            m = sum(n * p[n+1] for n in 0:nmax)
            v = sum((n - m)^2 * p[n+1] for n in 0:nmax)
            a_ins = sum(p[n+1] * min(1.0, lam / (n + 1)) for n in 0:nmax)
            a_del = sum(p[n+1] * min(1.0, n / lam) for n in 0:nmax) / (1 - p[1])
            return m, v, a_ins, a_del
        end

        function run_corner(z0V, seed; nburn=20_000, nsteps=200_000, thin=10)
            Random.seed!(seed)
            w = mkwalker()
            MC_grand_canonical_walk!(nburn, w, lj0, emax_inf; z0V=z0V, species=:Ar)
            ns = Int[]
            ins_att = 0; ins_acc = 0; del_att = 0; del_acc = 0
            for _ in 1:(nsteps ÷ thin)
                _, _, _, stats = MC_grand_canonical_walk!(thin, w, lj0, emax_inf; z0V=z0V, species=:Ar)
                push!(ns, w.list_num_par[1])
                ins_att += stats.insert_attempted; ins_acc += stats.insert_accepted
                del_att += stats.delete_attempted; del_acc += stats.delete_accepted
            end
            return mean(ns), var(ns), count(iszero, ns) / length(ns),
                   ins_acc / max(ins_att, 1), del_acc / max(del_att, 1), w
        end

        @testset "acceptance-ratio helpers (pre-move N convention)" begin
            gir = MonteCarloMoves.gc_insert_acceptance_ratio
            gdr = MonteCarloMoves.gc_delete_acceptance_ratio
            @test gir(8.0, 0, 0.25, 0.25) == 8.0
            @test gir(8.0, 7, 0.25, 0.25) == 1.0
            @test gir(8.0, 15, 0.25, 0.25) == 0.5
            @test gdr(8.0, 8, 0.25, 0.25) == 1.0
            @test gdr(8.0, 4, 0.25, 0.25) == 0.5
            @test gdr(8.0, 1, 0.25, 0.25) == 0.125
            @test gdr(8.0, 0, 0.25, 0.25) == 0.0
            # asymmetric channels carry the proposal-probability ratio
            @test gir(8.0, 3, 0.1, 0.4) == 8.0 * 4.0 / 4
            @test gdr(8.0, 3, 0.1, 0.4) == 0.25 * 3 / 8.0
        end

        @testset "argument validation" begin
            w = mkwalker()
            @test_throws ArgumentError MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=8.0, species=:Ar, p_move=0.9, p_insert=0.2)
            @test_throws ArgumentError MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=0.0, species=:Ar)
            wf = mkwalker()
            wf.frozen = [true]
            @test_throws ArgumentError MC_grand_canonical_walk!(1, wf, lj0, emax_inf; z0V=8.0, species=:Ar)
            tilted = [[10.0, 0.0, 0.0], [2.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
            wt = AtomWalker{1}(FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], tilted, pbc)))
            @test_throws ArgumentError MC_grand_canonical_walk!(1, wt, lj0, emax_inf; z0V=8.0, species=:Ar)
        end

        @testset "surgery helpers round-trip byte-identity" begin
            w = mkwalker()
            insert_particle!(w, SVector(1.0, 2.0, 3.0)u"Å", :Ar)
            insert_particle!(w, SVector(4.0, 5.0, 6.0)u"Å", :Ar)
            pos_before = copy(w.configuration.position)
            spc_before = copy(w.configuration.species)
            mss_before = copy(w.configuration.mass)
            insert_particle!(w, SVector(7.0, 8.0, 9.0)u"Å", :Ar)
            @test w.list_num_par == [3]
            @test length(w.configuration) == 3
            remove_particle!(w, 3)
            @test w.list_num_par == [2]
            @test w.configuration.position == pos_before
            @test w.configuration.species == spc_before
            @test w.configuration.mass == mss_before
            # order-preserving removal from the middle
            remove_particle!(w, 1)
            @test w.configuration.position == [pos_before[2]]
            @test_throws BoundsError remove_particle!(w, 5)
        end

        @testset "RNG-stream contract (forced branches)" begin
            # Insertion, ratio >= 1: accept with NO Metropolis draw (4 draws total)
            Random.seed!(778)
            probe = [rand() for _ in 1:8]
            @test probe[1] < 0.75    # branch precondition: channel draw lands in the insertion window
            Random.seed!(778)
            w = mkwalker()
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=1.0e6, species=:Ar, p_move=0.0, p_insert=0.75)
            @test w.list_num_par == [1]
            @test (stats.insert_attempted, stats.insert_accepted) == (1, 1)
            @test rand() == probe[5]

            # Insertion, ratio < 1: Metropolis uniform IS drawn (5 draws; p_delete = 0 forces ratio 0)
            Random.seed!(778)
            w = mkwalker()
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=8.0, species=:Ar, p_move=0.0, p_insert=1.0)
            @test w.list_num_par == [0]
            @test (stats.insert_attempted, stats.insert_accepted) == (1, 0)
            @test rand() == probe[6]

            # Deletion, ratio >= 1: accept with NO Metropolis draw (2 draws)
            Random.seed!(777)
            probe = [rand() for _ in 1:8]
            @test probe[1] >= 0.25   # branch precondition: channel draw lands in the deletion window
            Random.seed!(777)
            w = mkwalker()
            insert_particle!(w, SVector(5.0, 5.0, 5.0)u"Å", :Ar)
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=0.01, species=:Ar, p_move=0.0, p_insert=0.25)
            @test w.list_num_par == [0]
            @test (stats.delete_attempted, stats.delete_accepted) == (1, 1)
            @test rand() == probe[3]

            # Deletion, ratio < 1: Metropolis uniform IS drawn (3 draws)
            Random.seed!(777)
            w = mkwalker()
            insert_particle!(w, SVector(5.0, 5.0, 5.0)u"Å", :Ar)
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=1.0e6, species=:Ar, p_move=0.0, p_insert=0.25)
            @test w.list_num_par == [1]
            @test (stats.delete_attempted, stats.delete_accepted) == (1, 0)
            @test rand() == probe[4]

            # Ceiling rejection draws NO Metropolis uniform even at ratio < 1
            # (the ceiling is checked first): 4 draws total, walker reverted
            Random.seed!(778)
            probe778 = [rand() for _ in 1:8]
            Random.seed!(778)
            w = mkwalker()
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, 0.0u"eV"; z0V=0.3, species=:Ar, p_move=0.0, p_insert=0.75)
            @test w.list_num_par == [0]
            @test (stats.insert_attempted, stats.insert_accepted) == (1, 0)
            @test rand() == probe778[5]

            # Guard skips consume only the channel draw and are not counted as attempts
            Random.seed!(777)
            w = mkwalker()
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=8.0, species=:Ar, p_move=0.0, p_insert=0.25)
            @test w.list_num_par == [0]
            @test stats == (move_attempted=0, move_accepted=0, insert_attempted=0,
                            insert_accepted=0, delete_attempted=0, delete_accepted=0)
            @test rand() == probe[2]
            Random.seed!(779)
            probe779 = [rand() for _ in 1:4]
            Random.seed!(779)
            w = mkwalker()
            _, _, _, stats = MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=8.0, species=:Ar, p_move=1.0, p_insert=0.0)
            @test stats.move_attempted == 0
            @test rand() == probe779[2]
        end

        @testset "two-state birth-death closure at n_max = 1" begin
            Random.seed!(42424)
            w = mkwalker()
            MC_grand_canonical_walk!(10_000, w, lj0, emax_inf; z0V=0.6, species=:Ar, n_max=1)
            n1 = 0
            tot = 400_000
            for _ in 1:tot
                MC_grand_canonical_walk!(1, w, lj0, emax_inf; z0V=0.6, species=:Ar, n_max=1)
                n1 += w.list_num_par[1]
            end
            @test abs(n1 / tot - 0.6 / 1.6) < 0.013
            @test w.list_num_par[1] <= 1
        end

        @testset "ideal-gas stationary law (seeded, calibrated)" begin
            # (z0V, seed, tol_mean, tol_var, tol_p0, tol_ains, tol_adel)
            corners = [(0.5, 42421, 0.012, 0.015, 0.012, 0.010, 1.0e-12),
                       (8.0, 42422, 0.10, 1.0, 4.1e-4, 0.005, 0.027),
                       (32.0, 42423, 1.05, 3.4, Inf, 0.011, 0.018)]
            for (z0V, seed, tol_m, tol_v, tol_p0, tol_ai, tol_ad) in corners
                m_ref, v_ref, ai_ref, ad_ref = poisson_refs(z0V)
                m, v, p0, ai, ad, w = run_corner(z0V, seed)
                @test abs(m - m_ref) < tol_m
                @test abs(v - v_ref) < tol_v
                if isfinite(tol_p0)
                    @test abs(p0 - poisson_pmf(z0V, 0)) < tol_p0
                end
                @test abs(ai - ai_ref) < tol_ai
                @test abs(ad - ad_ref) < tol_ad
                # the configuration arrays and the particle count never drift apart
                @test length(w.configuration) == w.list_num_par[1]
                @test length(w.configuration.species) == w.list_num_par[1]
            end
            # sensitivity: the mean gate at the smallest corner detects the classic
            # off-by-one insertion bug, whose exact stationary mean follows from the
            # birth-death recursion p(n+1)/p(n) = a_ins(n)/a_del(n+1)
            lam = 0.5
            logp = zeros(101)
            for n in 0:99
                num = min(1.0, lam / max(n, 1))       # off-by-one foil
                den = min(1.0, (n + 1) / lam)
                logp[n+2] = logp[n+1] + log(num) - log(den)
            end
            pfoil = exp.(logp .- maximum(logp))
            pfoil ./= sum(pfoil)
            foil_shift = abs(sum(n * pfoil[n+1] for n in 0:100) - lam)
            @test 0.012 < 0.5 * foil_shift
        end

        @testset "fixed-seed digit pins (stream shape)" begin
            Random.seed!(424242)
            w = mkwalker()
            MC_grand_canonical_walk!(5000, w, lj0, emax_inf; z0V=8.0, species=:Ar)
            @test w.list_num_par[1] == 9
            @test w.energy == 0.0u"eV"
            lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5)
            Random.seed!(424243)
            w2 = mkwalker()
            MC_grand_canonical_walk!(5000, w2, lj, 1.0e5u"eV"; z0V=4.0, species=:Ar, step_size=1.0)
            @test w2.list_num_par[1] == 4
            # rtol 1e-12: tolerate compiler-level rounding drift along the same
            # trajectory (observed stable across architectures at authoring); a
            # stream or trajectory change fails loudly
            @test isapprox(ustrip(u"eV", w2.energy), -0.0014312826195886537; rtol=1e-12)
        end
    end

end


@testset "atomistic muVT kernel tests" begin
    using Random
    # Calibration ledger (gates at >= 3x the max three-seed deviation, per stat):
    # - Ideal-gas Poisson closure, zV = 3, T = 300 K, seeds 91001/91002/91003:
    #   max devs mean 0.036, var 0.059, p0 0.0018; gates 0.11, 0.18, 0.0055.
    # - Fixed-N pair-distance Boltzmann ratio (bins [2.5,2.9)/[3.3,3.7) A,
    #   quadrature reference 0.7929), seeds 92001/92002/92003: max rel dev
    #   0.0039; gate 0.012 relative.
    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    pbc = (true, true, true)
    seed_at = FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc))
    mkempty() = FastSystem(cell_vectors(seed_at), periodicity(seed_at),
                           empty(position(seed_at, :)), empty(species(seed_at, :)),
                           empty(mass(seed_at, :)))
    mkwalker() = (w = AtomWalker{1}(mkempty()); w.energy = 0.0u"eV"; w)
    lj0 = LJParameters(epsilon=0.0)
    T300 = 300.0
    logfact(n) = n == 0 ? 0.0 : sum(log, 1:n)
    poisson_pmf(lam, n) = exp(n * log(lam) - lam - logfact(n))

    @testset "argument validation" begin
        w = mkwalker()
        @test_throws ArgumentError MC_muVT_walk!(1, w, lj0, T300; zV=0.0, species=:Ar)
        @test_throws ArgumentError MC_muVT_walk!(1, w, lj0, 0.0; zV=1.0, species=:Ar)
        @test_throws ArgumentError MC_muVT_walk!(1, w, lj0, T300; zV=1.0, species=:Ar, p_move=0.7, p_insert=0.4)
        frozen_w = AtomWalker(FastSystem(atomic_system([:Ar => [1.0, 1.0, 1.0]u"Å"], box, pbc));
                              freeze_species=[:Ar])
        @test_throws ArgumentError MC_muVT_walk!(1, frozen_w, lj0, T300; zV=1.0, species=:Ar)
    end

    @testset "RNG-stream contract (forced branches)" begin
        # Insertion, combined ratio >= 1 (zV huge, zero-interaction): accept
        # with NO Metropolis draw (channel + three position uniforms = 4)
        Random.seed!(778)
        probe = [rand() for _ in 1:8]
        @test probe[1] < 0.75
        Random.seed!(778)
        w = mkwalker()
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=1.0e6, species=:Ar, p_move=0.0, p_insert=0.75)
        @test w.list_num_par == [1]
        @test (stats.insert_attempted, stats.insert_accepted) == (1, 1)
        @test rand() == probe[5]

        # Insertion, combined < 1: the Metropolis uniform IS drawn (5 draws;
        # p_delete = 0 forces the ratio to zero)
        Random.seed!(778)
        w = mkwalker()
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=8.0, species=:Ar, p_move=0.0, p_insert=1.0)
        @test w.list_num_par == [0]
        @test (stats.insert_attempted, stats.insert_accepted) == (1, 0)
        @test rand() == probe[6]

        # Deletion, combined >= 1: accept with NO Metropolis draw (2 draws)
        Random.seed!(777)
        probe = [rand() for _ in 1:8]
        @test probe[1] >= 0.25
        Random.seed!(777)
        w = mkwalker()
        insert_particle!(w, SVector(5.0, 5.0, 5.0)u"Å", :Ar)
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=0.01, species=:Ar, p_move=0.0, p_insert=0.25)
        @test w.list_num_par == [0]
        @test (stats.delete_attempted, stats.delete_accepted) == (1, 1)
        @test rand() == probe[3]

        # Deletion, combined < 1: the Metropolis uniform IS drawn (3 draws)
        Random.seed!(777)
        w = mkwalker()
        insert_particle!(w, SVector(5.0, 5.0, 5.0)u"Å", :Ar)
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=1.0e6, species=:Ar, p_move=0.0, p_insert=0.25)
        @test w.list_num_par == [1]
        @test (stats.delete_attempted, stats.delete_accepted) == (1, 0)
        @test rand() == probe[4]

        # Zero-interaction displacement: dU = 0 accepts downhill-style with NO
        # Metropolis draw (channel + index + three walk displacements = 5)
        Random.seed!(779)
        probe779 = [rand() for _ in 1:8]
        @test probe779[1] < 1.0
        Random.seed!(779)
        w = mkwalker()
        insert_particle!(w, SVector(5.0, 5.0, 5.0)u"Å", :Ar)
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=1.0, species=:Ar, p_move=1.0, p_insert=0.0)
        @test (stats.move_attempted, stats.move_accepted) == (1, 1)
        @test rand() == probe779[6]

        # Uphill displacement from the exact pair minimum: any move raises the
        # pair energy, so the Metropolis uniform IS drawn (6 draws)
        lj_up = LJParameters(epsilon=0.05, sigma=2.5, cutoff=3.0, shift=true)
        r0 = 2.5 * 2.0^(1 / 6)
        Random.seed!(781)
        probe781 = [rand() for _ in 1:8]
        Random.seed!(781)
        wp = AtomWalker(FastSystem(atomic_system(
            [:Ar => [2.0, 5.0, 5.0]u"Å", :Ar => [2.0 + r0, 5.0, 5.0]u"Å"], box, pbc)))
        wp.energy = interacting_energy(wp.configuration, lj_up, wp.list_num_par, wp.frozen)
        _, _, stats = MC_muVT_walk!(1, wp, lj_up, T300; zV=1.0, species=:Ar,
                                    p_move=1.0, p_insert=0.0, step_size=0.3)
        @test stats.move_attempted == 1
        @test rand() == probe781[7]

        # Guard skips consume only the channel draw and count no attempt
        Random.seed!(777)
        w = mkwalker()
        _, _, stats = MC_muVT_walk!(1, w, lj0, T300; zV=8.0, species=:Ar, p_move=0.0, p_insert=0.25)
        @test stats == (move_attempted=0, move_accepted=0, insert_attempted=0,
                        insert_accepted=0, delete_attempted=0, delete_accepted=0)
        @test rand() == probe[2]
    end

    @testset "ideal-gas Poisson closure (seeded, calibrated)" begin
        Random.seed!(91001)
        w = mkwalker()
        MC_muVT_walk!(20_000, w, lj0, T300; zV=3.0, species=:Ar)
        ns = Int[]
        for _ in 1:20_000
            MC_muVT_walk!(10, w, lj0, T300; zV=3.0, species=:Ar)
            push!(ns, w.list_num_par[1])
        end
        @test abs(mean(ns) - 3.0) < 0.11
        @test abs(var(ns) - 3.0) < 0.18
        @test abs(count(iszero, ns) / length(ns) - poisson_pmf(3.0, 0)) < 0.0055
        @test length(w.configuration) == w.list_num_par[1]
    end

    @testset "fixed-N Boltzmann law for the pair distance (seeded, calibrated)" begin
        lj = LJParameters(epsilon=0.02, sigma=2.5, cutoff=2.5, shift=true)
        kb = 8.617333262e-5
        β = 1 / (kb * T300)
        u_of(r) = ustrip(u"eV", FreeBird.AbstractPotentials.lj_energy((r)u"Å", lj))
        simpson(f, a, b; n=2000) = (h = (b - a) / n;
            h / 3 * sum(k -> f(a + k * h) * (k == 0 || k == n ? 1 : (isodd(k) ? 4 : 2)), 0:n))
        R_ref = simpson(r -> r^2 * exp(-β * u_of(r)), 2.5, 2.9) /
                simpson(r -> r^2 * exp(-β * u_of(r)), 3.3, 3.7)
        Random.seed!(92001)
        w = AtomWalker(FastSystem(atomic_system(
            [:Ar => [2.0, 5.0, 5.0]u"Å", :Ar => [4.8, 5.0, 5.0]u"Å"], box, pbc)))
        w.energy = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen)
        MC_muVT_walk!(50_000, w, lj, T300; zV=1.0, species=:Ar, p_move=1.0, p_insert=0.0, step_size=0.8)
        c1 = 0; c2 = 0
        for _ in 1:200_000
            MC_muVT_walk!(5, w, lj, T300; zV=1.0, species=:Ar, p_move=1.0, p_insert=0.0, step_size=0.8)
            d = ustrip(u"Å", FreeBird.EnergyEval.pbc_dist(position(w.configuration, 1),
                                                          position(w.configuration, 2),
                                                          w.configuration))
            (2.5 <= d < 2.9) && (c1 += 1)
            (3.3 <= d < 3.7) && (c2 += 1)
        end
        @test abs(c1 / c2 - R_ref) / R_ref < 0.012
    end

    @testset "incremental-energy drift audit" begin
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5, shift=true)
        Random.seed!(93001)
        w = mkwalker()
        MC_muVT_walk!(20_000, w, lj, T300; zV=6.0, species=:Ar, step_size=0.6)
        E_re = interacting_energy(w.configuration, lj, w.list_num_par, w.frozen) + w.energy_frozen_part
        @test abs(ustrip(u"eV", w.energy - E_re)) < 5e-10
        @test length(w.configuration) == w.list_num_par[1]
    end

    @testset "same-seed determinism" begin
        lj = LJParameters(epsilon=0.01, sigma=2.5, cutoff=2.5, shift=true)
        outs = map(1:2) do _
            Random.seed!(94001)
            w = mkwalker()
            _, rate, stats = MC_muVT_walk!(5_000, w, lj, T300; zV=4.0, species=:Ar)
            (w.list_num_par[1], ustrip(u"eV", w.energy), rate, stats)
        end
        @test outs[1] == outs[2]
    end
end
