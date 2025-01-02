@testset "Energy Evaluation Tests" begin
    @testset "atomistic_energies.jl tests" begin

        # Set up a basic periodic box (10Å × 10Å × 10Å)
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"

        @testset "pbc_dist function tests" begin

            # Two-particle system
            coor_list = [:H => [0.0, 0.0, 0.0], :H => [0.5, 0.5, 0.5]]
            two_particle_sys = periodic_system(coor_list, box, fractional=true)

            # Six-particle system
            coor_list1 = [:H => [0.2, 0.3, 0.5], :H => [0.8, 0.3, 0.5], 
                          :H => [0.5, 0.2, 0.3], :H => [0.7, 0.8, 0.3], 
                          :H => [0.5, 0.5, 0.1], :H => [0.4, 0.3, 0.1]]
            six_particle_sys = periodic_system(coor_list1, box, fractional=true)

            @testset "Regular Distance" begin
                # Test points within the box (no boundary crossing)
                pos1 = [2.0, 2.0, 2.0]u"Å"
                pos2 = [4.0, 4.0, 4.0]u"Å"
                dist = pbc_dist(pos1, pos2, two_particle_sys)
                expected_dist = sqrt(12.0)u"Å"
                @test dist ≈ expected_dist atol=1e-10u"Å"
            end

            @testset "Periodic Boundary Crossing" begin
                # Test distance across periodic boundary
                pos1 = [1.0, 1.0, 1.0]u"Å"
                pos2 = [9.0, 1.0, 1.0]u"Å"
                dist = pbc_dist(pos1, pos2, two_particle_sys)
                @test dist ≈ 2.0u"Å" atol=1e-10u"Å"
        
                # Test distance across multiple boundaries
                pos3 = [9.0, 9.0, 9.0]u"Å"
                pos4 = [1.0, 1.0, 1.0]u"Å"
                dist = pbc_dist(pos3, pos4, two_particle_sys)
                @test dist ≈ sqrt(12.0)u"Å" atol=1e-10u"Å"
            end

            @testset "Edge Cases" begin
                # Test zero distance
                pos_same = [5.0, 5.0, 5.0]u"Å"
                @test pbc_dist(pos_same, pos_same, two_particle_sys) ≈ 0.0u"Å" atol=1e-10u"Å"
        
                # Test maximum possible distance in periodic system
                pos1 = [0.0, 0.0, 0.0]u"Å"
                pos2 = [5.0, 5.0, 5.0]u"Å"
                max_dist = pbc_dist(pos1, pos2, two_particle_sys)
                @test max_dist ≈ sqrt(75.0)u"Å" atol=1e-10u"Å"
            end

            @testset "Multi-particle System Tests" begin
                # Get positions from sys1
                pos1 = position(six_particle_sys.particles[1])
                pos2 = position(six_particle_sys.particles[2])
                pos3 = position(six_particle_sys.particles[3])
                pos4 = position(six_particle_sys.particles[4])
                pos5 = position(six_particle_sys.particles[5])
                pos6 = position(six_particle_sys.particles[6])
        
                # Test periodic conditions
                @test pbc_dist(pos1, pos2, six_particle_sys) == 4.0u"Å"
                @test pbc_dist(pos3, pos4, six_particle_sys) ≈ 4.472135u"Å" rtol=1e-6
                @test pbc_dist(pos5, pos6, six_particle_sys) ≈ 2.236068u"Å" rtol=1e-6
        
                # Test non-periodic conditions
                boundary_conditions = [DirichletZero(), DirichletZero(), DirichletZero()]
                non_periodic_sys = FlexibleSystem(six_particle_sys; boundary_conditions=boundary_conditions)
        
                @test pbc_dist(pos1, pos2, non_periodic_sys) ≈ 6.0u"Å"
                @test pbc_dist(pos3, pos4, non_periodic_sys) ≈ 6.324555u"Å" rtol=1e-6
                @test pbc_dist(pos5, pos6, non_periodic_sys) ≈ 2.236068u"Å" rtol=1e-6
            end

        end


        @testset "Energy evaluation functions tests" begin

            # System 1: 6 particles
            coor_list1 = [:H => [0.2, 0.3, 0.5], :H => [0.8, 0.3, 0.5], 
                          :H => [0.5, 0.2, 0.3], :H => [0.7, 0.8, 0.3],
                          :H => [0.5, 0.5, 0.1], :H => [0.4, 0.3, 0.1]]
            sys1 = periodic_system(coor_list1, box, fractional=true)

            # System 2: 3 particles
            coor_list2 = [:H => [0, 0, 0.55], :H => [0.3, 0.6, 0.6], :H => [0.5, 0.5, 0.65]]
            sys2 = periodic_system(coor_list2, box, fractional=true)

            # System 3: 6 particles in pairs
            coor_list3 = [:H => [0.2, 0.2, 0.3], :H => [0.4, 0.4, 0.3],
                          :H => [0.4, 0.4, 0.6], :H => [0.6, 0.6, 0.6],
                          :H => [0.6, 0.6, 0.9], :H => [0.8, 0.8, 0.9]]
            sys3 = periodic_system(coor_list3, box, fractional=true)
    
            # Setup LJ parameters
            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
            ljs1 = [LJParameters(epsilon=e) for e in [11, 21, 31, 12, 22, 32, 13, 23, 33]]
            ljs2 = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]
            ljp1 = CompositeLJParameters(3, ljs1)
            ljp2 = CompositeLJParameters(3, ljs2)

            @testset "inter_component_energy function tests" begin
                # Test basic energy calculation
                @test EnergyEval.inter_component_energy(sys1, sys2, lj) ≈ -0.2524365u"eV" rtol=1e-6
                
                # Test symmetry
                energy_12 = EnergyEval.inter_component_energy(sys1, sys2, lj)
                energy_21 = EnergyEval.inter_component_energy(sys2, sys1, lj)
                @test energy_12 ≈ energy_21 rtol=1e-10
                
                # Test empty system
                empty_sys = periodic_system([], box, fractional=true)
                @test EnergyEval.inter_component_energy(empty_sys, sys1, lj) == 0.0u"eV"
            end

            @testset "intra_component_energy function tests" begin
                @test EnergyEval.intra_component_energy(sys1, lj) ≈ 0.617334u"eV" rtol=1e-6
                @test EnergyEval.intra_component_energy(sys2, lj) ≈ 0.457052u"eV" rtol=1e-6
                
                # Test single particle system (should be zero)
                single_particle = periodic_system([:H => [0.0, 0.0, 0.0]], box, fractional=true)
                @test EnergyEval.intra_component_energy(single_particle, lj) == 0.0u"eV"
                
                # Test empty system
                empty_sys = periodic_system([], box, fractional=true)
                @test EnergyEval.intra_component_energy(empty_sys, lj) == 0.0u"eV"
            end

            @testset "Single LJ frozen_energy function tests" begin
                # Test with various frozen configurations
                @test frozen_energy(sys3, lj, [2, 2, 2], [true, false, false]) ≈ -0.099785393u"eV" rtol=1e-7
                @test frozen_energy(sys3, lj, [2, 2, 2], [true, true, false]) ≈ -0.32785795u"eV" rtol=1e-7
                
                # Test all frozen
                @test frozen_energy(sys3, lj, [2, 2, 2], [true, true, true]) ≈ -0.56551325u"eV" rtol=1e-7
                
                # Test none frozen
                @test frozen_energy(sys3, lj, [2, 2, 2], [false, false, false]) == 0.0u"eV"
                
                # Test error handling
                @test_throws ArgumentError frozen_energy(sys3, lj, [2, 2], [true, false, false])
            end

            @testset "Composite LJ frozen_energy function tests" begin
                @test frozen_energy(sys3, ljp1, [2, 2, 2], [true, true, false]) ≈ -0.49376980u"eV" rtol=1e-7
                @test frozen_energy(sys3, ljp2, [2, 2, 2], [true, true, false]) ≈ -0.42928048u"eV" rtol=1e-7
                
                # Test component mismatch error
                @test_throws ArgumentError frozen_energy(sys3, ljp1, [2, 2], [true, false])
            end

            @testset "interacting_energy functions tests" begin
                # Test single LJ potential
                @test interacting_energy(sys3, lj, [2, 2, 2], [true, false, false]) ≈ -0.46572786u"eV" rtol=1e-7
                
                # Test composite LJ potential
                @test interacting_energy(sys3, ljp1, [2, 2, 2], [true, false, false]) ≈ -0.68481194u"eV" rtol=1e-7
                
                # Test consistency with intra_component_energy for single component
                @test interacting_energy(sys3, lj) == EnergyEval.intra_component_energy(sys3, lj)
                @test interacting_energy(sys3, lj) ≈ -0.56551325u"eV" rtol=1e-7
                
                # Test empty system
                empty_sys = periodic_system([], box, fractional=true)
                @test interacting_energy(empty_sys, lj) == 0.0u"eV"
            end
        end
    end


    @testset "lattice_energies.jl tests" begin
        # Setup test parameters
        n_sites = 4
        neighbors = [
            [[2,4], [3]], # site 1's first and second neighbors
            [[1,3], [4]], # site 2's first and second neighbors
            [[2,4], [1]], # site 3's first and second neighbors
            [[1,3], [2]]  # site 4's first and second neighbors
        ]
        
        ham = GenericLatticeHamiltonian(-0.04u"eV", [-0.01, -0.0025]*u"eV")
        hams = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
        mlham = MLatticeHamiltonian(2, hams)

        # Setup lattice systems
        lattice = MLattice{1,SquareLattice}(
            lattice_constant=1.0,
            supercell_dimensions=(6,6,1),
            periodicity = (true, true, false),
            cutoff_radii = [1.1, 1.5],
            components = [[1, 2, 4]],
            adsorptions = [1, 2, 3]
        )

        @testset "lattice_interaction_energy function tests" begin
            empty_occupations = fill(false, n_sites)
            @test EnergyEval.lattice_interaction_energy(empty_occupations, neighbors, ham) ≈ 0.0u"eV"
    
            single_site = [true, false, false, false]
            @test EnergyEval.lattice_interaction_energy(single_site, neighbors, ham) ≈ 0.0u"eV"
    
            two_nn_sites = [true, true, false, false]
            @test EnergyEval.lattice_interaction_energy(two_nn_sites, neighbors, ham) ≈ ham.nth_neighbor_interactions[1]
    
            two_nnn_sites = [true, false, true, false]
            @test EnergyEval.lattice_interaction_energy(two_nnn_sites, neighbors, ham) ≈ ham.nth_neighbor_interactions[2]
    
            full_lattice = fill(true, n_sites)
            @test EnergyEval.lattice_interaction_energy(full_lattice, neighbors, ham) ≈ 
                (2*ham.nth_neighbor_interactions[1] + ham.nth_neighbor_interactions[2]) * 4 / 2
        end
    
        @testset "inter_component_energy function tests" begin
            lattice1 = [true, false, false, false]
            lattice2 = [false, true, false, false]
            @test EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham) ≈ ham.nth_neighbor_interactions[1]
    
            lattice1 = [true, false, true, false]
            lattice2 = [false, true, false, true]
            @test EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham) ≈ 4 * ham.nth_neighbor_interactions[1]
    
            # Test symmetry
            lattice1 = [true, false, false, false]
            lattice2 = [true, true, true, true]
            @test EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham) == 
                  EnergyEval.inter_component_energy(lattice2, lattice1, neighbors, ham)
        end
    
        @testset "interacting_energy functions tests" begin
            # Test empty lattice
            empty_lattice = deepcopy(lattice)
            empty_lattice.components[1] .= false
            @test interacting_energy(empty_lattice, ham) ≈ 0.0u"eV"
            @test interacting_energy(empty_lattice, mlham) ≈ 0.0u"eV"
    
            # Test single site
            single_site = deepcopy(lattice)
            single_site.components[1] .= false
            single_site.components[1][1] = true
            @test interacting_energy(single_site, ham) ≈ ham.on_site_interaction
            @test interacting_energy(single_site, mlham) ≈ ham.on_site_interaction
    
            # Test adjacent sites
            adjacent_sites = deepcopy(lattice)
            adjacent_sites.components[1] .= false
            adjacent_sites.components[1][1:2] .= true
            @test interacting_energy(adjacent_sites, ham) ≈ 
                  2 * ham.on_site_interaction + ham.nth_neighbor_interactions[1]
            @test interacting_energy(adjacent_sites, mlham) ≈ 
                  2 * ham.on_site_interaction + ham.nth_neighbor_interactions[1]
        end
    end


    @testset "helpers.jl tests" begin
        # Setup test system
        box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
        
        # Mixed species system
        coor_list = [:H => [0.1, 0.1, 0.1], :H => [0.2, 0.2, 0.2],  # H atoms
                     :O => [0.3, 0.3, 0.3],                          # O atom
                     :Fe => [0.4, 0.4, 0.4], :Fe => [0.5, 0.5, 0.5], # Fe atoms
                     :Au => [0.6, 0.6, 0.6], :Au => [0.7, 0.7, 0.7]  # Au atoms
                    ]

        mixed_system = periodic_system(coor_list, box, fractional=true)

        # Single species system for empty component tests
        single_particle = periodic_system([:H => [0.0, 0.0, 0.0]], box, fractional=true)

        @testset "split_components function tests" begin
            components = split_components(mixed_system, [2, 1, 2, 2])
            @test length(components) == 4
            @test length(components[1]) == 2  # H component
            @test length(components[2]) == 1  # O component
            @test length(components[3]) == 2  # Fe component
            @test length(components[4]) == 2  # Au component
            
            # Test single particle system
            single_comp = split_components(single_particle, [1])
            @test length(single_comp) == 1
            @test length(single_comp[1]) == 1
        end

        @testset "split_components_by_chemical_species function" begin
            components = EnergyEval.split_components_by_chemical_species(mixed_system)
            @test length(components) == 4  # H, O, Fe, Au
            
            # Verify component sizes
            species_counts = [length(comp) for comp in components]
            @test sort(species_counts) == [1, 2, 2, 2]
            
            # Test single particle system
            single_comps = EnergyEval.split_components_by_chemical_species(single_particle)
            @test length(single_comps) == 1
            @test length(single_comps[1]) == 1
        end

        @testset "check_num_components" begin
            # Valid case
            @test nothing === EnergyEval.check_num_components(3, [2, 3, 1], [true, false, true])
            
            # Error cases
            @test_throws ArgumentError EnergyEval.check_num_components(3, [2, 3], [true, false, true])
            @test_throws ArgumentError EnergyEval.check_num_components(3, [2, 3, 1], [true, false])
        end

        @testset "sort_components_by_atomic_number" begin
            # Test with merge_same_species=true
            list_num_par, new_sys = sort_components_by_atomic_number(mixed_system, merge_same_species=true)
            @test length(list_num_par) == 4
            @test sum(list_num_par) == length(mixed_system)
            @test length(new_sys) == length(mixed_system)
            
            # Test with merge_same_species=false
            list_num_par, new_sys = sort_components_by_atomic_number(mixed_system, merge_same_species=false)
            @test length(list_num_par) == 4
            @test sum(list_num_par) == length(mixed_system)
            
            # Test single particle system
            list_num_par, new_sys = sort_components_by_atomic_number(single_particle)
            @test length(list_num_par) == 1
            @test length(new_sys) == 1
        end
    end    

end