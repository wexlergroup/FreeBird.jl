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

@testset "partition function and internal energy functions tests" begin
    
    @testset "Basic functionality" begin
        # Error cases of incorrect applying
        @test_throws MethodError partition_function(1.0, 1.0, 0.0)
        @test_throws MethodError partition_function(0, [1.0], [0.0])
        @test_throws MethodError partition_function(1, [-1.0], [0.0])
        @test_throws MethodError partition_function(1, [1.0], [-1.0])
        @test_throws DimensionMismatch partition_function(1.0, [0.5, 0.5], [0.0, 1.0, 10.0])
        @test_throws DimensionMismatch partition_function(1.0, [0.3, 0.3, 0.4], [0.0, 1.0])

        @test_throws MethodError internal_energy(1.0, 1.0, 0.0)
        @test_throws MethodError internal_energy(0, [1.0], [0.0])
        @test_throws MethodError internal_energy(1, [-1.0], [0.0])
        @test_throws MethodError internal_energy(1, [1.0], [-1.0])
        @test_throws DimensionMismatch internal_energy(1.0, [0.5, 0.5], [0.0, 1.0, 10.0])
        @test_throws DimensionMismatch internal_energy(1.0, [0.3, 0.3, 0.4], [0.0, 1.0])


        # Correct single-term case
        @test partition_function(1.0, [1.0], [0.0]) ≈ 1.0
        @test internal_energy(1.0, [1.0], [0.0]) == 0.0


        # multiple-term cases
        ωi = collect(0.40:-0.1:0.10)
        Ei = collect(0.0:0.3:1.0)
        @test partition_function(1.0, ωi, Ei) ≈ 0.7726647593973806
        @test internal_energy(1.0, ωi, Ei) ≈ 0.2188818676047736

    end

    @testset "Temperature dependence" begin
        # Test high temperature limit (β ➡ 0)
        @test partition_function(1e-10, ωi, Ei) ≈ 0.99999999997
        @test internal_energy(1e-10, ωi, Ei) ≈ 0.299999999991


        # Test low temperature limit (β ➡ ∞)
        β_large = 1.0e10
        ground_energy_idx = argmin(Ei)
        @test partition_function(β_large, ωi, Ei) ≈ ωi[ground_energy_idx] * exp(-Ei[ground_energy_idx]*β_large)
        @test internal_energy(β_large, ωi, Ei) ≈ ωi[ground_energy_idx] * Ei[ground_energy_idx] * exp(-Ei[ground_energy_idx]*β_large) / ωi[ground_energy_idx] * exp(-Ei[ground_energy_idx]*β_large)

    end

    @testset "Edge cases " begin
        # Empty system
        @test partition_function(1.0, Float64[], Float64[]) == 0.0
        @test isnan(internal_energy(1.0, Float64[], Float64[]))
        
        # Zero weights
        @test partition_function(1.0, [0.0], [1.0]) == 0.0
        @test isnan(internal_energy(1.0, [0.0], [1.0]))
        
        # Very large energies (potential overflow)
        @test partition_function(1.0, [1.0], [1e3]) ≈ 0.0
        @test isnan(internal_energy(1.0, [1.0], [1e3]))
        
        # Very small energies (potential underflow)
        @test isinf(partition_function(1.0, [1.0], [-1e3]))
        @test isnan(internal_energy(1.0, [1.0], [-1e3]))

    end

    @testset "Numerical stability" begin
        # Test with very large β
        β_large = 1e20
        weights = [0.9, 0.1]
        energies = [0.0, 1.0]
        result_Z = partition_function(β_large, weights, energies)
        result_U = internal_energy(β_large, weights, energies)
        @test !isnan(result_Z)
        @test !isinf(result_Z)
        @test !isnan(result_U)
        @test !isinf(result_U)

        
        # Test with very large energies
        energies_large = [1e4, 1e4 + 1.0]
        result_Z = partition_function(1.0, weights, energies_large)
        result_U = internal_energy(1.0, weights, energies_large)
        @test !isnan(result_Z)
        @test !isinf(result_Z)
        @test isnan(result_U)

        
        # Test preservation of relative scales
        β = 1.0
        weights = [0.5, 0.5]
        E1 = [0.0, 1e-10]
        E2 = [0.0, 1e-9]
        Z1 = partition_function(β, weights, E1)
        Z2 = partition_function(β, weights, E2)
        @test Z1 ≈ Z2 rtol=1e-8
    end

end
