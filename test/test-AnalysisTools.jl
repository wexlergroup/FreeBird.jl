@testset "AnalysisTools.jl tests" begin
    @testset "ωᵢ function tests" begin
        # Basic tests
        @test ωᵢ([1], 4) ≈ [1/5 * (4/5)^1]
        @test ωᵢ([1, 2, 3], 4) ≈ [1/5 * (4/5)^1, 1/5 * (4/5)^2, 1/5 * (4/5)^3]
        @test ωᵢ([0], 4) ≈ [1/5]
        @test ωᵢ([100], 4)[1] ≈ 1/5 * (4/5)^100 rtol=1e-10
        
        # Edge cases and properties
        @test ωᵢ(Int[], 4) == Float64[]
        @test_throws MethodError ωᵢ([1.5], 4)
        result = ωᵢ(collect(1:5), 4)
        @test all(diff(result) .< 0)  # Monotonic decrease
        @test all(result .> 0)        # All positive
    end

    @testset "Partition function and internal energy tests" begin
        # Basic functionality
        @test partition_function(1.0, [1.0], [0.0]) ≈ 1.0
        @test internal_energy(1.0, [1.0], [0.0]) == 0.0

        # Temperature limits
        ωi = collect(0.40:-0.1:0.10)
        Ei = collect(0.0:0.3:1.0)
        @test partition_function(1e-10, ωi, Ei) ≈ 1.0 rtol=1e-8
        @test internal_energy(1e-10, ωi, Ei) ≈ 0.3 rtol=1e-8

        # Edge cases
        @test partition_function(1.0, Float64[], Float64[]) == 0.0
        @test isnan(internal_energy(1.0, Float64[], Float64[]))
        @test_throws DimensionMismatch partition_function(1.0, [0.5, 0.5], Ei)
    end

    @testset "Constant-volume heat capacity tests" begin
        kb = 8.617333262e-5
        ωi = [0.4, 0.3, 0.2, 0.1]
        Ei = [0.0, 0.3, 0.6, 0.9]
        dof = 3

        # Basic functionality and limits
        @test cv(1.0, [1.0], [0.0], 3) ≈ 3 * kb / 2.0 rtol=1e-8
        @test cv(1e-10, ωi, Ei, dof) ≈ dof * kb / 2.0 rtol=1e-8
        
        # DataFrame interface
        df = DataFrame(iter = [1, 2, 3], emax = [1.0, 1.5, 2.0])
        β = 1.0
        βs = [0.1, 1.0, 10.0]
        n_walkers = 4
        result = cv(df, βs, dof, n_walkers)
        @test length(result) == length(βs)
        @test all(.!isnan.(result))

        # Edge cases
        @test cv(1.0, ωi, Ei, 0) ≥ 0
        @test isnan(cv(β, Float64[], Float64[], dof))
        @test_throws DimensionMismatch cv(1.0, [0.5, 0.5], Ei, dof)
    end
end