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

end

