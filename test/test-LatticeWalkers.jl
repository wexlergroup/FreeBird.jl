@testset "Lattice Walkers" begin
    # 4x4 square lattice
    dims = (2, 2, 1)
    sl = SLattice{SquareLattice}(supercell_dimensions=dims)
    @test sl.neighbors == [[[2, 3], [4]], [[1, 4], [3]], [[1, 4], [2]], [[2, 3], [1]]]
    @test sl.components == [[true, true, true, true]]
    @test sl.basis == [(0.0, 0.0, 0.0)]
    @test sl.positions == [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 1.0 0.0]
end