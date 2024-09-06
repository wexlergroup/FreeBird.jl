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

    @test EnergyEval.inter_component_energy(sys1, sys2, lj) ≈ -0.5382944918043817u"eV"

    
    @test EnergyEval.intra_component_energy(sys1, lj) ≈ -0.2597893542152926u"eV"
    @test EnergyEval.intra_component_energy(sys2, lj) ≈ -0.00023134752229080925u"eV"



end