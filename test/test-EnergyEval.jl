# @testset "Energy Evaluation Tests" begin

#     @testset "EnergyEval.jl" begin
        
#         # Generate different periodic systems
#         box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"

#         coor_list1 = [:H => [0.2, 0.3, 0.5], :H => [0.8, 0.3, 0.5], :H => [0.5, 0.2, 0.3], :H => [0.7, 0.8, 0.3], :H => [0.5, 0.5, 0.1], :H => [0.4, 0.3, 0.1]]
#         sys1 = periodic_system(coor_list1, box, fractional=true)

#         coor_list2 = [:H => [0, 0, 0.55], :H => [0.3, 0.6, 0.6], :H => [0.5, 0.5, 0.65]]
#         sys2 = periodic_system(coor_list2, box, fractional=true)
#         coor_list3 = [:H => [0.2, 0.2, 0.3], :H => [0.4, 0.4, 0.3], :H => [0.4, 0.4, 0.6], :H => [0.6, 0.6, 0.6], :H => [0.6, 0.6, 0.9], :H => [0.8, 0.8, 0.9]]
#         sys3 = periodic_system(coor_list3, box, fractional=true)

#         # Generate single LJ potential.
#         lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

#         # Generate composite LJ potential
#         ljs1 = [LJParameters(epsilon=e) for e in [11, 21, 31, 12, 22, 32, 13, 23, 33]]
#         ljs2 = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]

#         ljp1 = CompositeLJParameters(3, ljs1)
#         ljp2 = CompositeLJParameters(3, ljs2)

        
#         @testset "pbc_dist function" begin

#             pos1 = position(sys1.particles[1])
#             pos2 = position(sys1.particles[2])
#             pos3 = position(sys1.particles[3])
#             pos4 = position(sys1.particles[4])
#             pos5 = position(sys1.particles[5])
#             pos6 = position(sys1.particles[6])

#             # Test output distances in x, y, z directions under the periodic condition
#             @test pbc_dist(pos1, pos2, sys1) == 4.0u"Å"
#             @test pbc_dist(pos3, pos4, sys1) ≈ 4.472135u"Å"           rtol=1e-6
#             @test pbc_dist(pos5, pos6, sys1) ≈ 2.236068u"Å"           rtol=1e-6


#             # Test output distances in x, y, z directions under the non-periodic condition
#             boundary_conditions = [DirichletZero(), DirichletZero(), DirichletZero()]
#             non_periodic_sys = FlexibleSystem(sys1; boundary_conditions=boundary_conditions)

#             @test pbc_dist(pos1, pos2, non_periodic_sys) ≈ 6.0u"Å"
#             @test pbc_dist(pos3, pos4, non_periodic_sys) ≈ 6.324555u"Å"         rtol=1e-6
#             @test pbc_dist(pos5, pos6, sys1) ≈ 2.236068u"Å"                     rtol=1e-6

#         end

#         @testset "Energy evaluation functions" begin

#             # Test inter_component_energy() function for calculating the energy between two components of a system
#             @test EnergyEval.inter_component_energy(sys1, sys2, lj) ≈ -0.2524365u"eV"        rtol=1e-6

#             # Test threaded implementation against non-threaded
#             EnergyEval.inter_component_energy(sys1, sys2, lj) ≈ EnergyEval.inter_component_energy_threaded(sys1, sys2, lj)

#             # Test intra_component_energy() function for calculating the energy within a component of a system
#             @test EnergyEval.intra_component_energy(sys1, lj) ≈ 0.617334u"eV"               rtol=1e-6
#             @test EnergyEval.intra_component_energy(sys2, lj) ≈ 0.457052u"eV"               rtol=1e-6
            
#             # Test frozen_energy() function for calculating the frozen energy in a system using the single LJ potential.
#             @test frozen_energy(sys3, lj, [2, 2, 2], [true, false, false]) ≈ -0.099785393u"eV"       rtol=1e-7

#             # Test frozen_energy() function for calculating the frozen energy in a system using the composite LJ potential.
#             @test frozen_energy(sys3, ljp1, [2, 2, 2], [true, true, false]) ≈ -0.49376980u"eV"      rtol=1e-7
#             @test frozen_energy(sys3, ljp2, [2, 2, 2], [true, true, false]) ≈ -0.42928048u"eV"      rtol=1e-7

#             # Test whether interacting_energy() function is consistent with intra_component_energy() function when the single componenet system is not frozen
#             @test interacting_energy(sys3, lj) == EnergyEval.intra_component_energy(sys3, lj)

#             # Test interacting_energy() function for calculating the interacting energy in a system using the single LJ potential.
#             @test interacting_energy(sys3, lj, [2, 2, 2], [true, false, false]) ≈ -0.46572786u"eV"      rtol=1e-7

#             # Test interacting_energy() function for calculating the interacting energy in a system using the composite LJ potential.
#             @test interacting_energy(sys3, ljp1, [2, 2, 2], [true, false, false]) ≈ -0.68481194u"eV"    rtol=1e-7
#             @test interacting_energy(sys3, ljp2, [2, 2, 2], [true, false, false]) ≈ interacting_energy(sys3, ljp1, [2, 2, 2], [true, false, false])    # Shouldn't this test be false?

#             # Test for interacting_energy(lattice::LatticeSystem, h::LatticeGasHamiltonian) function
#             @test interacting_energy(sys3, lj) == EnergyEval.intra_component_energy(sys3, lj)
#             @test interacting_energy(sys3, lj) ≈ -0.56551325u"eV"            rtol=1e-7

#             empty_sys = periodic_system([], box, fractional=true)
#             @test interacting_energy(empty_sys, lj) == 0.0u"eV"

#         end

#         @testset "Energy evaluation functions for lattice models" begin

#             # Generate different Hamiltonians
#             ham = GenericLatticeHamiltonian(-0.04u"eV", [-0.01, -0.0025]*u"eV")
#             hams = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
#             mlham = MLatticeHamiltonian(2, hams)

#             # Generate lattice
#             lattice = MLattice{1,SquareLattice}(
#                 lattice_constant=1.0,
#                 supercell_dimensions=(6,6,1),
#                 periodicity = (true, true, false),
#                 cutoff_radii = [1.1, 1.5],
#                 components = [[1, 2, 4]],
#                 adsorptions = [1, 2, 3] 
#             )
            
#             # Setup square lattice with periodic boundaries
#             n_sites = 4
#             neighbors = [
#                 [[2,4], [3]], # site 1's first and second neighbors
#                 [[1,3], [4]], # site 2's first and second neighbors
#                 [[2,4], [1]], # site 3's first and second neighbors
#                 [[1,3], [2]]  # site 4's first and second neighbors
#             ]
           

#             # Test for interacting_energy(lattice::SLattice, h::LatticeGasHamiltonian) function
#             empty_lattice = deepcopy(lattice)
#             empty_lattice.components[1] .= false
#             @test EnergyEval.interacting_energy(empty_lattice, ham) ≈ 0.0u"eV"
#             @test EnergyEval.interacting_energy(empty_lattice, mlham) ≈ 0.0u"eV"

#             single_site = deepcopy(lattice)
#             single_site.components[1] .= false
#             single_site.components[1][1] = true
#             @test EnergyEval.interacting_energy(single_site, ham) ≈ ham.on_site_interaction
#             @test EnergyEval.interacting_energy(single_site, mlham) ≈ ham.on_site_interaction

#             adjacent_sites = deepcopy(lattice)
#             adjacent_sites.components[1] .= false
#             adjacent_sites.components[1][1:2] .= true
#             @test EnergyEval.interacting_energy(adjacent_sites, ham) ≈ 2 * ham.on_site_interaction + ham.nth_neighbor_interactions[1]
#             @test EnergyEval.interacting_energy(adjacent_sites, mlham) ≈ 2 * ham.on_site_interaction + ham.nth_neighbor_interactions[1]

#             # Test for lattice_interaction_energy() function
#             empty_occupations = fill(false, n_sites)
#             @test EnergyEval.lattice_interaction_energy(empty_occupations, neighbors, ham) ≈ 0.0u"eV"

#             single_site = [true, false, false, false]
#             @test EnergyEval.lattice_interaction_energy(single_site, neighbors, ham) ≈ 0.0u"eV"

#             two_nn_sites = [true, true, false, false]
#             @test EnergyEval.lattice_interaction_energy(two_nn_sites, neighbors, ham) ≈ ham.nth_neighbor_interactions[1]
        
#             two_nnn_sites = [true, false, true, false]
#             @test EnergyEval.lattice_interaction_energy(two_nnn_sites, neighbors, ham) ≈ ham.nth_neighbor_interactions[2]

#             full_lattice = fill(true, n_sites)
#             @test EnergyEval.lattice_interaction_energy(full_lattice, neighbors, ham) ≈ (2*ham.nth_neighbor_interactions[1] + 1*ham.nth_neighbor_interactions[2]) * 4 / 2

#             # Test for inter_component_energy() function
#             lattice1 = [true, false, false, false]
#             lattice2 = [false, true, false, false]
#             @test EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham) ≈ ham.nth_neighbor_interactions[1]

#             lattice1 = [true, false, true, false]
#             lattice2 = [false, true, false, true]
#             @test EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham) ≈ 4 * ham.nth_neighbor_interactions[1]

#             lattice1 = [true, false, false, false]
#             lattice2 = [true, true, true, true]
#             result1 = EnergyEval.inter_component_energy(lattice1, lattice2, neighbors, ham)
#             result2 = EnergyEval.inter_component_energy(lattice2, lattice1, neighbors, ham)
#             @test result1 == result2

#             # Test for interacting_energy() function

#         end

#     end

#     @testset "helpers.jl" begin
        
#     end

# end



@testset "Energy Evaluation Tests" begin
    @testset "atomistic_energies.jl tests" begin
        @testset "pbc_dist function tests" begin
            # Set up a basic periodic box (10Å × 10Å × 10Å)
            box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
            
            # Create a simple periodic system with two particles
            coor_list = [:H => [0.0, 0.0, 0.0], :H => [0.5, 0.5, 0.5]]
            sys = periodic_system(coor_list, box, fractional=true)

            @testset "Regular Distance" begin
                # Test points within the box (no boundary crossing)
                pos1 = [2.0, 2.0, 2.0]u"Å"
                pos2 = [4.0, 4.0, 4.0]u"Å"
                dist = pbc_dist(pos1, pos2, sys)
                expected_dist = sqrt(12.0)u"Å"
                @test dist ≈ expected_dist atol=1e-10u"Å"
            end

            @testset "Periodic Boundary Crossing" begin
                # Test distance across periodic boundary
                pos1 = [1.0, 1.0, 1.0]u"Å"
                pos2 = [9.0, 1.0, 1.0]u"Å"
                dist = pbc_dist(pos1, pos2, sys)
                @test dist ≈ 2.0u"Å" atol=1e-10u"Å"
        
                # Test distance across multiple boundaries
                pos3 = [9.0, 9.0, 9.0]u"Å"
                pos4 = [1.0, 1.0, 1.0]u"Å"
                dist = pbc_dist(pos3, pos4, sys)
                @test dist ≈ sqrt(12.0)u"Å" atol=1e-10u"Å"
            end

            @testset "Edge Cases" begin
                # Test zero distance
                pos_same = [5.0, 5.0, 5.0]u"Å"
                @test pbc_dist(pos_same, pos_same, sys) ≈ 0.0u"Å" atol=1e-10u"Å"
        
                # Test maximum possible distance in periodic system
                pos1 = [0.0, 0.0, 0.0]u"Å"
                pos2 = [5.0, 5.0, 5.0]u"Å"
                max_dist = pbc_dist(pos1, pos2, sys)
                @test max_dist ≈ sqrt(75.0)u"Å" atol=1e-10u"Å"
            end
        end


        @testset "Energy evaluation functions tests" begin

            # Set up basic test systems
            box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    
            # Create simple Lennard-Jones parameters
            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)

            @testset "Intra-Component Energy Tests" begin
                # Test system with two particles
                sys_two = periodic_system(
                    [:H => [0.0, 0.0, 0.0]u"Å", :H => [2.5, 0.0, 0.0]u"Å"],
                    box
                )
                
                # Test single particle system (should have zero energy)
                sys_single = periodic_system(
                    [:H => [0.0, 0.0, 0.0]u"Å"],
                    box
                )
        
                # Calculate energies
                e_two = EnergyEval.intra_component_energy(sys_two, lj)
                e_single = EnergyEval.intra_component_energy(sys_single, lj)
        
                # Tests
                @test typeof(e_two) == typeof(0.0u"eV")
                @test e_single == 0.0u"eV"  # Single particle should have zero energy
                @test e_two < 0.0u"eV"      # Particles at sigma distance should have negative energy
            end

        end

    end

end