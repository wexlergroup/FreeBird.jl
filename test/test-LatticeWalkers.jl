@testset "Lattice Walkers" begin
    # test the neighbors function
    dims = (4, 4)

    # 4x4 square lattice
    @test sort(AbstractWalkers.neighbors(:square, dims, 1, 1)) == [(1, 2), (1, 4), (2, 1), (4, 1)]  # corner site
    @test sort(AbstractWalkers.neighbors(:square, dims, 1, 2)) == [(1, 1), (1, 3), (2, 2), (4, 2)]  # edge site
    @test sort(AbstractWalkers.neighbors(:square, dims, 2, 2)) == [(1, 2), (2, 1), (2, 3), (3, 2)]  # center site

    # 4x4 hexagonal lattice

    # corner sites
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 1, 1)) == [(1, 2), (1, 4), (2, 1), (2, 2), (4, 1), (4, 2)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 1, 4)) == [(1, 1), (1, 3), (2, 1), (2, 4), (4, 1), (4, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 4, 1)) == [(1, 1), (1, 4), (3, 1), (3, 4), (4, 2), (4, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 4, 4)) == [(1, 3), (1, 4), (3, 3), (3, 4), (4, 1), (4, 3)]

    # edge sites
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 1, 2)) == [(1, 1), (1, 3), (2, 2), (2, 3), (4, 2), (4, 3)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 2, 1)) == [(1, 1), (1, 4), (2, 2), (2, 4), (3, 1), (3, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 1, 3)) == [(1, 2), (1, 4), (2, 3), (2, 4), (4, 3), (4, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 3, 1)) == [(2, 1), (2, 2), (3, 2), (3, 4), (4, 1), (4, 2)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 2, 4)) == [(1, 3), (1, 4), (2, 1), (2, 3), (3, 3), (3, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 4, 2)) == [(1, 1), (1, 2), (3, 1), (3, 2), (4, 1), (4, 3)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 3, 4)) == [(2, 1), (2, 4), (3, 1), (3, 3), (4, 1), (4, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 4, 3)) == [(1, 2), (1, 3), (3, 2), (3, 3), (4, 2), (4, 4)]

    # center sites
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 2, 2)) == [(1, 1), (1, 2), (2, 1), (2, 3), (3, 1), (3, 2)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 3, 3)) == [(2, 3), (2, 4), (3, 2), (3, 4), (4, 3), (4, 4)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 2, 3)) == [(1, 2), (1, 3), (2, 2), (2, 4), (3, 2), (3, 3)]
    @test sort(AbstractWalkers.neighbors(:hexagonal, dims, 3, 2)) == [(2, 2), (2, 3), (3, 1), (3, 3), (4, 2), (4, 3)]
end