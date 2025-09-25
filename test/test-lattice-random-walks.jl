@testset "Monte Carlo Moves tests" begin
    @testset "lattice random walk tests" begin

        sl = SLattice{SquareLattice}(components=[[1,2,3,4]])
        new_sl = lattice_random_walk!(deepcopy(sl))

        @test occupied_site_count(sl) == occupied_site_count(new_sl)
        @test length(sl.components[1]) == length(new_sl.components[1])

        ml = MLattice{2,SquareLattice}()
        new_ml = lattice_random_walk!(deepcopy(ml))

        @test occupied_site_count(ml) == occupied_site_count(new_ml)
        @test length(ml.components[1]) == length(new_ml.components[1])
        @test length(ml.components[2]) == length(new_ml.components[2])
    end
end
