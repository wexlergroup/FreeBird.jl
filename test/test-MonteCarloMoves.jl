using Test
using StaticArrays
using Random
using Unitful

@testset "Monte Carlo Moves tests" begin
    
    @testset "random_walks.jl tests" begin

        @testset "single_atom_random_walk function" begin

            # Test basic functionality
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 2.0, 3.0] .* u"Å")
            step_size = 0.5
            new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
            
            @test new_pos ≠ pos
            
            # Test step size constraints
            for _ in 1:10
                pos = SVector{3, typeof(0.0u"Å")}([0.0, 0.0, 0.0] .* u"Å")
                new_pos = MonteCarloMoves.single_atom_random_walk!(pos, step_size)
                
                # Check that no component moved more than step_size
                @test all(abs.(new_pos - pos) .<= step_size * u"Å")
            end
            
            # Test type stability
            pos = SVector{3, typeof(0.0u"Å")}([1.0, 1.0, 1.0] .* u"Å")
            @test typeof(MonteCarloMoves.single_atom_random_walk!(pos, step_size)) == typeof(pos)

        end

        @testset "MC_random_walk functions" begin

            @testset "AtomWalker MC_random_walk" begin
                # Setup test parameters
                n_steps = 100
                step_size = 0.5
                emax = 1.0u"eV"
                
                # Create a simple test system
                config = Configuration() # Create a test configuration
                at = AtomWalker(config, 0.0u"eV", [1], [false])  # Example values
                lj = LJParameters()  # Create test LJ parameters
                
                # Run MC simulation
                accept_this_walker, accept_rate, updated_at = MC_random_walk!(
                    n_steps, at, lj, step_size, emax
                )
                
                # Tests
                @test typeof(accept_this_walker) == Bool
                @test 0.0 <= accept_rate <= 1.0
                @test typeof(updated_at) == typeof(at)
                @test updated_at.energy <= emax
            end

            @testset "LatticeWalker MC_random_walk" begin

            end




        end
        

    end

    @testset "nve_walks.jl tests" begin

    end

    @testset "helpers.jl tests" begin

    end

end

