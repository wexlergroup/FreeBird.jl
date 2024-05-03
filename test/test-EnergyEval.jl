@testset "Energy Evaluation" begin
    list = [:H => [0.2, 0.5, 0.5],:H => [0.8, 0.5, 0.5]]
    box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
    sys = periodic_system(list,box,fractional=true)
    pos1 = position(sys.particles[1])
    pos2 = position(sys.particles[2])
    @test pbc_dist(pos1, pos2, sys) == 4.0u"Å"
    boundary_conditions = [DirichletZero(), DirichletZero(), DirichletZero()]
    sys = FlexibleSystem(sys; boundary_conditions=boundary_conditions)
    @test pbc_dist(pos1, pos2, sys) == 6.0u"Å"
end