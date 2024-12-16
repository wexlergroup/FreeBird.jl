using Test

@testset "ωᵢ function tests" begin
    @testset "Basic functionality" begin
        # Test with single iteration
        @test ωᵢ([1], 4) ≈ [1/5 * (4/5)^1]
        
        # Test with multiple iterations
        result = ωᵢ([1, 2, 3], 4)
        @test length(result) == 3
        @test result ≈ [1/5 * (4/5)^1, 1/5 * (4/5)^2, 1/5 * (4/5)^3]
    end

    @testset "Edge cases" begin
        # Test with zero iterations
        @test ωᵢ([0], 4) ≈ [1/5]
        
        # Test with large iterations
        large_iter = 100
        result = ωᵢ([large_iter], 4)
        @test result[1] ≈ 1/5 * (4/5)^large_iter rtol=1e-10
        
        # Test with single walker
        result = ωᵢ([1, 2, 3], 1)
        @test result ≈ [1/2 * (1/2)^1, 1/2 * (1/2)^2, 1/2 * (1/2)^3]
    end

    @testset "Mathematical properties" begin
        # Test monotonic decrease
        result = ωᵢ(collect(1:5), 4)
        for i in 1:length(result)-1
            @test result[i] > result[i+1]
        end
        
        # Test all values are positive
        @test all(result .> 0)
        
        # Test with different n_walkers values
        n_walkers_list = [1, 2, 5, 10, 100]
        for n in n_walkers_list
            result = ωᵢ([1], n)
            @test result[1] ≈ (1/(n+1))*(n/(n+1))^1
        end
    end

    @testset "Vector properties" begin
        # Test empty vector
        @test ωᵢ(Int[], 4) == Float64[]
        
        # Test non-integer iterations
        @test_throws MethodError ωᵢ([1.5], 4)
        
        # Test negative iterations
        result = ωᵢ([-1], 4)
        @test result[1] ≈ 1/5 * (4/5)^(-1)
    end

    @testset "Numerical stability" begin
        # Test with very large n_walkers
        large_n = 10^6
        result = ωᵢ([1], large_n)
        @test !isnan(result[1])
        @test result[1] > 0
        
        # Test sum of weights for consecutive iterations
        n_walkers = 4
        iters = collect(1:10)
        weights = ωᵢ(iters, n_walkers)
        total_weight = sum(weights)
        @test total_weight < 1  # Sum should be less than 1
        @test total_weight ≈ (1/(n_walkers+1)) * (n_walkers/(n_walkers+1))^iters[1] * (1-(n_walkers/(n_walkers+1))^iters[end])/(1-(n_walkers/(n_walkers+1)))
    end
end