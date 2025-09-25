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

    @testset "Swap occupied sites across components tests" begin
        ml = MLattice{2,SquareLattice}(components=[[1,3],[2,4]])
        new_ml = deepcopy(ml)
        MonteCarloMoves.swap_occupied_sites_across_components!(new_ml, 1, 2)

        @test new_ml.components[1][1] == 0
        @test new_ml.components[1][2] == 1
        @test new_ml.components[2][1] == 1
        @test new_ml.components[2][2] == 0
    end

    @testset "Swap empty and occupied sites tests" begin
        ml = MLattice{2,SquareLattice}(components=[[1],[2]])
        new_ml = deepcopy(ml)
        MonteCarloMoves.swap_empty_occupied_sites!(new_ml, 1, 2)

        @test new_ml.components[1][1] == 0
        @test new_ml.components[1][2] == 1
        @test new_ml.components[2][1] == 1
        @test new_ml.components[2][2] == 0
    end
end
