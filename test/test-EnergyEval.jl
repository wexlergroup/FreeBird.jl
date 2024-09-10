@testset "Energy Evaluation" begin

    coor_list1 = [:H => [0.2, 0.5, 0.5], :H => [0.8, 0.5, 0.5], :H => [0.5, 0.2, 0.5], :H => [0.5, 0.8, 0.5], :H => [0.5, 0.5, 0.2], :H => [0.5, 0.5, 0.8]]
    coor_list2 = [:H => [0, 0, 0], :H => [0.5, 0.5, 0.5]]
    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"

    sys1 = periodic_system(coor_list1, box, fractional=true)
    pos1 = position(sys1.particles[1])
    pos2 = position(sys1.particles[2])

    pos3 = position(sys1.particles[3])
    pos4 = position(sys1.particles[4])
    
    pos5 = position(sys1.particles[5])
    pos6 = position(sys1.particles[6])

    # Test output distances in x, y, z directions under the periodic condition
    @test pbc_dist(pos1, pos2, sys1) == 4.0u"Å" 
    @test pbc_dist(pos3, pos4, sys1) == 4.0u"Å"
    @test pbc_dist(pos5, pos6, sys1) == 4.0u"Å"

    boundary_conditions = [DirichletZero(), DirichletZero(), DirichletZero()]
    non_periodic_sys = FlexibleSystem(sys1; boundary_conditions=boundary_conditions)

    # Test output distances in x, y, z directions under the non-periodic condition
    @test pbc_dist(pos1, pos2, non_periodic_sys) == 6.0u"Å"
    @test pbc_dist(pos3, pos4, non_periodic_sys) == 6.0u"Å"
    @test pbc_dist(pos5, pos6, non_periodic_sys) == 6.0u"Å"

    
    sys2 = periodic_system(coor_list2, box, fractional=true)
    lj = LJParameters(epsilon=0.1,sigma=2.5,cutoff=3.5,shift=false)

    # Test inter_component_energy() function for calculating the energy between two components of a system
    @test EnergyEval.inter_component_energy(sys1, sys2, lj) ≈ -0.5382944918043817u"eV"

    # Test intra_component_energy() function for calculating the energy within a component of a system
    @test EnergyEval.intra_component_energy(sys1, lj) ≈ -0.2597893542152926u"eV"
    @test EnergyEval.intra_component_energy(sys2, lj) ≈ -0.00023134752229080925u"eV"


    coor_list3 = [:H => [0.2, 0.2, 0.3], :H => [0.4, 0.4, 0.3], :H => [0.4, 0.4, 0.6], :H => [0.6, 0.6, 0.6], :H => [0.6, 0.6, 0.9], :H => [0.8, 0.8, 0.9]]
    sys3 = periodic_system(coor_list3, box, fractional=true)
    
    @test EnergyEval.frozen_energy(sys3, lj, [2, 2, 2], [true, false, false]) ≈ -0.09978539310395718u"eV"

    ljs1 = [LJParameters(epsilon=e) for e in [11, 21, 31, 12, 22, 32, 13, 23, 33]]
    ljs2 = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]

    ljp1 = CompositeLJParameters(3, ljs1)
    ljp2 = CompositeLJParameters(3, ljs2)

    @test EnergyEval.frozen_energy(sys3, ljp1, [2, 2, 2], [true, true, false]) ≈ -0.49376979751528466u"eV"
    @test EnergyEval.frozen_energy(sys3, ljp2, [2, 2, 2], [true, true, false]) ≈ -0.4292804781822809u"eV"


end