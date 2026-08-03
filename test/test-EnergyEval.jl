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
                boundary_conditions = (false, false, false)
                non_periodic_sys = FastSystem(six_particle_sys.particles, box, boundary_conditions)
        
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
            ljp1 = CompositeParameterSets(3, ljs1)
            ljp2 = CompositeParameterSets(3, ljs2)

            @testset "pair energy function tests" begin
                
                dist = pbc_dist(position(sys1.particles[1]), position(sys1.particles[2]), sys1)
                pair_energy = EnergyEval.pair_energy(dist, lj)
                @test pair_energy ≈ -0.02242077243863605u"eV" rtol=1e-6
            end

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

                # Test error handling
                @test_throws ArgumentError interacting_energy(sys3, ljp1, [2, 2], [true, false, false])
            end

            @testset "single site energy function tests" begin
                # Test single LJ potential
                @test single_site_energy(1, sys3, lj, [length(sys3)]) ≈ -0.12184883176851288u"eV" rtol=1e-7

                # Test composite LJ potential
                @test single_site_energy(1, sys3, ljp1, [length(sys3)]) ≈ -0.09615777513293565u"eV" rtol=1e-7
                @test single_site_energy(2, sys3, ljp1, [3,3]) ≈ -0.1597175167721998u"eV" rtol=1e-7


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

        @testset "triangular lattice energies" begin
            using Random
            # First triangular coverage on the energy path: six first-shell
            # neighbors per site under PBC, and energy = J × (occupied
            # first-shell pairs) against an independent pair counter.
            J = 1.0
            ham_nn = GenericLatticeHamiltonian(0.0, [J], u"eV")
            tri = MLattice{1,TriangularLattice}(
                lattice_constant=1.0,
                supercell_dimensions=(3, 3, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1],
                components=:equal,
                adsorptions=:full)
            M = length(tri.components[1])
            @test M == 18
            @test all(length(tri.neighbors[i][1]) == 6 for i in 1:M)
            # First-shell adjacency is symmetric
            @test all(i in tri.neighbors[j][1] for i in 1:M for j in tri.neighbors[i][1])

            pair_count(occ) = sum(
                occ[i] && occ[j] for i in 1:M for j in tri.neighbors[i][1] if j > i)

            work = deepcopy(tri)
            work.components[1] .= false
            @test interacting_energy(work, ham_nn) ≈ 0.0u"eV"
            work.components[1][1] = true
            @test interacting_energy(work, ham_nn) ≈ 0.0u"eV"
            work.components[1][first(tri.neighbors[1][1])] = true
            @test interacting_energy(work, ham_nn) ≈ J * u"eV"
            work.components[1] .= true
            @test interacting_energy(work, ham_nn) ≈ J * 6 * M / 2 * u"eV"

            Random.seed!(42)
            for _ in 1:50
                occ = rand(Bool, M)
                work.components[1] .= occ
                @test interacting_energy(work, ham_nn) ≈ J * pair_count(occ) * u"eV"
            end

            # The customary two-shell cutoffs [1.1, 1.5] leave the second
            # shell EMPTY on a triangular lattice (second-neighbor distance
            # √3 ≈ 1.732) — the constructor warns about the silent shell.
            tri15 = @test_logs (:warn, r"empty on every") match_mode=:any MLattice{1,TriangularLattice}(
                lattice_constant=1.0,
                supercell_dimensions=(4, 4, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.5],
                components=:equal,
                adsorptions=:full)
            @test all(isempty(tri15.neighbors[i][2])
                      for i in 1:length(tri15.components[1]))
            # ... while [1.1, 1.8] captures the six √3 second neighbors.
            tri18 = MLattice{1,TriangularLattice}(
                lattice_constant=1.0,
                supercell_dimensions=(4, 4, 1),
                periodicity=(true, true, false),
                cutoff_radii=[1.1, 1.8],
                components=:equal,
                adsorptions=:full)
            M18 = length(tri18.components[1])
            @test all(length(tri18.neighbors[i][1]) == 6 for i in 1:M18)
            @test all(length(tri18.neighbors[i][2]) == 6 for i in 1:M18)
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

    @testset "multi-shell pair couplings (shell-count validation, sign-mixed sets)" begin
        using Random

        function shell_lattice(L, cutoffs; image_multiplicity::Bool=false)
            MLattice{1,SquareLattice}(
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(L, L, 1),
                periodicity=(true, true, false),
                cutoff_radii=cutoffs,
                components=[[false for _ in 1:L*L]],
                adsorptions=:full,
                image_multiplicity=image_multiplicity)
        end

        # Square-lattice shells at distances 1, √2, 2, √5, 2√2, 3, √10, √13, 4
        cut4 = [1.1, 1.5, 2.1, 2.3]
        cut9 = [1.1, 1.5, 2.1, 2.3, 2.9, 3.1, 3.2, 3.7, 4.1]

        # ------------------------------------------------------------
        @testset "Zhang pair part: closed-form adlayer energies (8x8, 4 shells)" begin
            # Pair part of the Zhang, Blum & Reuter O-Pd(100) lattice-gas
            # Hamiltonian [PRB 75, 235406 (2007)], FreeBird sign convention
            # (positive = repulsive); four shells with mixed signs
            V0 = -1.249
            J = [0.292, 0.090, -0.050, -0.010]
            ham = GenericLatticeHamiltonian(V0, J, u"eV")
            lat = shell_lattice(8, cut4)

            # c(2x2) checkerboard at half coverage: each occupied site sees 4
            # occupied 2nd-shell (√2) and 4 occupied 3rd-shell (distance 2)
            # neighbors; shells 1 and 4 connect opposite sublattices and
            # contribute exactly zero
            lat.components[1] .= [iseven(((s - 1) % 8) + ((s - 1) ÷ 8)) for s in 1:64]
            @test interacting_energy(lat, ham).val ≈
                  32 * (V0 + 2 * J[2] + 2 * J[3]) atol = 1e-10

            # p(2x2) at quarter coverage: only the 3rd shell (distance 2)
            # connects occupied sites
            lat.components[1] .= [((s - 1) % 8) % 2 == 0 && ((s - 1) ÷ 8) % 2 == 0
                                  for s in 1:64]
            @test interacting_energy(lat, ham).val ≈
                  16 * (V0 + 2 * J[3]) atol = 1e-10
        end

        # ------------------------------------------------------------
        @testset "Kappus-shape nine-shell set vs independent brute force (16x16)" begin
            # Kappus [Surf. Sci. 691, 121444 (2020)] elastic O-Pd(100) pair
            # set: nine shells, mixed signs. The 16x16 cell is faithful for
            # every shell (circumference 16 > 2 x 4).
            Random.seed!(138)
            J9 = [0.228, 0.066, -0.043, -0.015, 0.020, 0.007, -0.006, -0.009, 0.008]
            ham9 = GenericLatticeHamiltonian(0.0, J9, u"eV")
            L = 16
            lat = shell_lattice(L, cut9)

            # Independent reference: integer minimum-image squared distances
            # classify the shells exactly; shares no code with
            # compute_neighbors. (r² = 17, i.e. √17 ≈ 4.12, falls outside the
            # 4.1 cutoff, matching the library's shell assignment.)
            shell_of_r2 = Dict(1 => 1, 2 => 2, 4 => 3, 5 => 4, 8 => 5,
                               9 => 6, 10 => 7, 13 => 8, 16 => 9)
            function brute_force_E(occ)
                E = 0.0
                for a in 1:L*L, b in (a+1):L*L
                    (occ[a] && occ[b]) || continue
                    ia, ja = (a - 1) % L, (a - 1) ÷ L
                    ib, jb = (b - 1) % L, (b - 1) ÷ L
                    dx = min(mod(ia - ib, L), mod(ib - ia, L))
                    dy = min(mod(ja - jb, L), mod(jb - ja, L))
                    sh = get(shell_of_r2, dx^2 + dy^2, 0)
                    sh == 0 || (E += J9[sh])
                end
                return E
            end

            # atol sized for summation-order headroom (~2000 terms, |E| up
            # to ~28 eV: observed discrepancy ~4e-13); adjacent shell
            # couplings differ by >= 1e-3, so 1e-9 loses no bug-catching power
            for _ in 1:50
                occ = rand(L * L) .< 0.4
                lat.components[1] .= occ
                @test interacting_energy(lat, ham9).val ≈ brute_force_E(occ) atol = 1e-9
            end
        end

        # ------------------------------------------------------------
        @testset "shell-count mismatch behavior" begin
            lat4 = shell_lattice(4, cut4)
            lat4.components[1][1:2] .= true

            # More coupled shells than the lattice provides: ArgumentError,
            # not a raw BoundsError — on every energy path and at liveset
            # construction (where energies are assigned)
            ham5 = GenericLatticeHamiltonian(0.0, [0.1, 0.1, 0.1, 0.1, 0.1], u"eV")
            @test_throws ArgumentError interacting_energy(lat4, ham5)
            @test_throws ArgumentError FreeBird.EnergyEval.lattice_interaction_energy(
                lat4.components[1], lat4.neighbors, ham5)
            @test_throws ArgumentError FreeBird.EnergyEval.inter_component_energy(
                lat4.components[1], lat4.components[1], lat4.neighbors, ham5)
            walkers5 = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
            @test_throws ArgumentError LatticeGasWalkers(walkers5, ham5)

            # Fewer coupled shells than the lattice carries: legal (outer
            # shells simply uncoupled), but warned once at liveset construction
            ham2 = GenericLatticeHamiltonian(0.0, [0.1, 0.05], u"eV")
            E2 = interacting_energy(lat4, ham2)
            @test unit(E2) == u"eV"
            walkers2 = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
            @test_logs (:warn, r"couples only the first 2") match_mode = :any LatticeGasWalkers(
                walkers2, ham2)

            # Matching counts: silent
            ham4 = GenericLatticeHamiltonian(0.0, [0.1, 0.05, 0.02, 0.01], u"eV")
            walkers4 = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
            @test_logs min_level = Base.CoreLogging.Warn LatticeGasWalkers(walkers4, ham4)

            # Both guards fire through the MLatticeHamiltonian dispatch too
            # (its shell count is the shared type parameter N, not C)
            mham2 = MLatticeHamiltonian(2,
                [GenericLatticeHamiltonian(0.0, [0.1, 0.05], u"eV") for _ in 1:3])
            lat1 = shell_lattice(4, [1.1])
            lat1.components[1][1:2] .= true
            @test_throws ArgumentError interacting_energy(lat1, mham2)
            walkers_m = [LatticeWalker(deepcopy(lat4), energy=0.0u"eV", iter=0)]
            @test_logs (:warn, r"couples only the first 2") match_mode = :any LatticeGasWalkers(
                walkers_m, mham2)

            # The warn also reaches the raw-lattice sampler entry points,
            # which never construct a liveset
            wl_params = WangLandauParameters(energy_min=-0.5, energy_max=3.0,
                                             num_energy_bins=30, num_steps=10,
                                             max_iter=2, random_seed=7)
            @test_logs (:warn, r"couples only the first 2") match_mode = :any wang_landau(
                deepcopy(lat4), ham2, wl_params)
            @test_logs (:warn, r"couples only the first 2") match_mode = :any nvt_monte_carlo(
                MCNewSample(), deepcopy(lat4), ham2, 300.0, 5, 7)
        end

        # ------------------------------------------------------------
        @testset "sign-mixed multi-shell GC-NS smoke run" begin
            Random.seed!(1380)
            ham = GenericLatticeHamiltonian(-0.05, [0.292, 0.090, -0.050, -0.010], u"eV")
            lat = shell_lattice(8, cut4)
            walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                       for _ in 1:20]
            ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
            params = IdealGasReferencedGCNSParameters(
                mc_steps=20, reference_fugacity=1.0, energy_perturbation=1e-9)
            save = SaveEveryN("t_shell.csv", "t_shell.traj", "t_shell.ls",
                              1000000, 1000000, 1000000)
            df, final_ls, _ = ideal_gas_referenced_nested_sampling(
                ls, params, Int64(100),
                MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3), save)
            rm.(["t_shell.csv", "t_shell.traj", "t_shell.ls"], force=true)

            @test nrow(df) > 0
            # Live-set energies match direct recomputation up to the
            # tie-breaking perturbation
            for w in final_ls.walkers
                @test isapprox(w.energy.val,
                               interacting_energy(w.configuration, ham).val;
                               atol=1e-8)
            end
        end

        # ------------------------------------------------------------
        @testset "image multiplicity: faithful cells unchanged" begin
            # On cells faithful for every shell (circumference > 2 x outermost
            # cutoff) the two conventions must agree element-for-element, and
            # neither construction may warn
            for (L, cutoffs) in ((4, [1.1, 1.5]), (8, cut4), (16, cut9))
                lat_def = @test_logs min_level = Base.CoreLogging.Warn shell_lattice(L, cutoffs)
                lat_mult = @test_logs min_level = Base.CoreLogging.Warn shell_lattice(
                    L, cutoffs; image_multiplicity=true)
                @test lat_mult.neighbors == lat_def.neighbors
            end
        end

        # ------------------------------------------------------------
        @testset "image multiplicity: closed-form energies on wrapped cells" begin
            # 4x4 with a unit third-shell coupling: the third shell (distance
            # 2 = L/2) is halved under the minimum-image convention (16 eV),
            # while multiplicity recovers the bulk-tiled 32 eV
            ham3 = GenericLatticeHamiltonian(0.0, [0.0, 0.0, 1.0], u"eV")
            sq_def = @test_logs (:warn, r"wrap the periodic cell") match_mode = :any shell_lattice(
                4, [1.1, 1.5, 2.1])
            sq_def.components[1] .= true
            @test all(length.(sq_def.neighbors[i]) == [4, 4, 2] for i in 1:16)
            @test interacting_energy(sq_def, ham3).val ≈ 16.0 atol = 1e-10
            sq_mult = @test_logs min_level = Base.CoreLogging.Warn shell_lattice(
                4, [1.1, 1.5, 2.1]; image_multiplicity=true)
            sq_mult.components[1] .= true
            @test all(length.(sq_mult.neighbors[i]) == [4, 4, 4] for i in 1:16)
            @test interacting_energy(sq_mult, ham3).val ≈ 32.0 atol = 1e-10

            # Shipped triangular (4, 2, 1) cell with a unit second-shell
            # coupling: the 5-bond vs 6-bond per-site closed forms
            ham2 = GenericLatticeHamiltonian(0.0, [0.0, 1.0], u"eV")
            tri_def = @test_logs (:warn, r"wrap the periodic cell") match_mode = :any MLattice{1,TriangularLattice}()
            tri_def.components[1] .= true
            @test interacting_energy(tri_def, ham2).val ≈ 40.0 atol = 1e-10
            tri_mult = @test_logs min_level = Base.CoreLogging.Warn MLattice{1,TriangularLattice}(
                image_multiplicity=true)
            tri_mult.components[1] .= true
            @test interacting_energy(tri_mult, ham2).val ≈ 48.0 atol = 1e-10

            # Bulk-tiled chain closed forms, periodicity (true, false, false)
            chain(d1, cutoffs; kwargs...) = MLattice{1,SquareLattice}(;
                lattice_constant=1.0,
                basis=[(0.0, 0.0, 0.0)],
                supercell_dimensions=(d1, 1, 1),
                periodicity=(true, false, false),
                cutoff_radii=cutoffs,
                components=[[true for _ in 1:d1]],
                adsorptions=:full,
                kwargs...)

            # (4,1,1) with a unit second-shell coupling: one wrapped
            # second-shell bond per site by default, two when tiled
            hamc = GenericLatticeHamiltonian(0.0, [0.0, 1.0], u"eV")
            c4_def = @test_logs (:warn, r"wrap the periodic cell") match_mode = :any chain(4, [1.1, 2.1])
            @test interacting_energy(c4_def, hamc).val ≈ 2.0 atol = 1e-10
            c4_mult = chain(4, [1.1, 2.1]; image_multiplicity=true)
            @test interacting_energy(c4_mult, hamc).val ≈ 4.0 atol = 1e-10

            # (2,1,1) with J1 = J2 = 1: multiplicity counts the doubled
            # nearest-neighbor bond (2 eV, shell 1) plus the self-image bonds
            # (2 eV, shell 2) — the bulk-tiled chain value; the minimum-image
            # convention sees a single bond in total
            hamcc = GenericLatticeHamiltonian(0.0, [1.0, 1.0], u"eV")
            c2_def = @test_logs (:warn, r"wrap the periodic cell") match_mode = :any chain(2, [1.1, 2.1])
            @test interacting_energy(c2_def, hamcc).val ≈ 1.0 atol = 1e-10
            c2_mult = chain(2, [1.1, 2.1]; image_multiplicity=true)
            @test interacting_energy(c2_mult, hamcc).val ≈ 4.0 atol = 1e-10
        end

        # ------------------------------------------------------------
        @testset "multiplicity lattice through GC-NS end-to-end" begin
            Random.seed!(2081)
            ham = GenericLatticeHamiltonian(-0.05, [0.292, 0.090, -0.050], u"eV")
            lat = shell_lattice(4, [1.1, 1.5, 2.1]; image_multiplicity=true)
            walkers = [LatticeWalker(deepcopy(lat), energy=0.0u"eV", iter=0)
                       for _ in 1:20]
            ls = LatticeGasWalkers(walkers, ham; assign_energy=false)
            params = IdealGasReferencedGCNSParameters(
                mc_steps=20, reference_fugacity=1.0, energy_perturbation=1e-9)
            save = SaveEveryN("t_mult.csv", "t_mult.traj", "t_mult.ls",
                              1000000, 1000000, 1000000)
            df, final_ls, _ = ideal_gas_referenced_nested_sampling(
                ls, params, Int64(100),
                MCGrandCanonicalMoves(p_move=0.4, p_insert=0.3), save)
            rm.(["t_mult.csv", "t_mult.traj", "t_mult.ls"], force=true)

            @test nrow(df) > 0
            # Live-set energies stay consistent with direct recomputation on
            # the duplicated neighbor lists, up to the tie-breaking
            # perturbation
            for w in final_ls.walkers
                @test isapprox(w.energy.val,
                               interacting_energy(w.configuration, ham).val;
                               atol=1e-8)
            end
        end
    end

end