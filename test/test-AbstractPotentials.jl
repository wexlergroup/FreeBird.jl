@testset "AbstractPotentials Tests" begin
    @testset "LJParameters struct and functions tests" begin
        @testset "Default constructor tests" begin
            lj = LJParameters()
            @test lj.epsilon == 1.0u"eV"
            @test lj.sigma == 1.0u"Å"
            @test lj.cutoff == Inf
            @test lj.shift == 0.0u"eV"
        end
    
        @testset "Constructor with custorm parameters tests" begin
            # Test with epsilon only
            lj = LJParameters(epsilon=0.5)
            @test lj.epsilon == 0.5u"eV"
            @test lj.sigma == 1.0u"Å"
            @test lj.cutoff == Inf
            @test lj.shift == 0.0u"eV"
    
            # Test with sigma only
            lj = LJParameters(sigma=2.5)
            @test lj.epsilon == 1.0u"eV"
            @test lj.sigma == 2.5u"Å"
            @test lj.cutoff == Inf
            @test lj.shift == 0.0u"eV"
    
            # Test with all parameters
            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
            @test lj.epsilon == 0.1u"eV"
            @test lj.sigma == 2.5u"Å"
            @test lj.cutoff == 3.5
            @test lj.shift == 0.0u"eV"
        end
    
        @testset "Shift parameter behavior tests" begin
            # Test with shift=true and finite cutoff
            lj = LJParameters(cutoff=3.5, shift=true)
            expected_shift = lj_energy(1.0u"eV", 1.0u"Å", 3.5u"Å")
            @test lj.shift ≈ expected_shift
    
            # Test with shift=false
            lj = LJParameters(cutoff=3.5, shift=false)
            @test lj.shift == 0.0u"eV"
    
            # Test with explicit shift value
            lj = LJParameters(cutoff=3.5, shift=5.0)
            @test lj.shift == 5.0u"eV"
    
            # Test with infinite cutoff
            lj = LJParameters(cutoff=Inf, shift=true)
            @test lj.shift == 0.0u"eV"
        end
    
        @testset "Physical parameter ranges tests" begin
            # Test with zero epsilon
            lj = LJParameters(epsilon=0.0)
            @test lj.epsilon == 0.0u"eV"
    
            # Test with very small sigma
            lj = LJParameters(sigma=1e-10)
            @test lj.sigma == 1e-10u"Å"
    
            # Test with very large cutoff
            lj = LJParameters(cutoff=1e10)
            @test lj.cutoff == 1e10
    
            # Test with negative shift
            lj = LJParameters(shift=-1.0)
            @test lj.shift == -1.0u"eV"
        end
    
        @testset "Unit consistency tests" begin
            # Test that units are properly assigned
            lj = LJParameters(epsilon=2.0, sigma=3.0)
            @test typeof(lj.epsilon) == typeof(1.0u"eV")
            @test typeof(lj.sigma) == typeof(1.0u"Å")
            @test typeof(lj.shift) == typeof(1.0u"eV")
            @test typeof(lj.cutoff) == Float64
        end
    
        @testset "Documentation examples tests" begin
            # Test all examples from the docstring
            lj = LJParameters(epsilon=0.1, sigma=2.5, cutoff=3.5, shift=false)
            @test lj.epsilon == 0.1u"eV"
            @test lj.sigma == 2.5u"Å"
            @test lj.cutoff == 3.5
            @test lj.shift == 0.0u"eV"
    
            lj = LJParameters(sigma=2.5)
            @test lj.sigma == 2.5u"Å"
            @test lj.cutoff == Inf
    
            lj = LJParameters(cutoff=3.5, shift=5.0)
            @test lj.cutoff == 3.5
            @test lj.shift == 5.0u"eV"
    
            lj = LJParameters(cutoff=3.5, shift=true)
            @test lj.cutoff == 3.5
            @test lj.shift ≈ -0.0021747803916549904u"eV" rtol=1e-10
    
            lj = LJParameters(cutoff=3.5, shift=false)
            @test lj.cutoff == 3.5
            @test lj.shift == 0.0u"eV"
        end
    
        @testset "Error handling tests" begin
            # Test handling of unspecified units (should be automatically added)
            lj = LJParameters(epsilon=2.0)
            @test lj.epsilon == 2.0u"eV"
            @test lj.sigma == 1.0u"Å"
    
            # Test type stability
            lj = LJParameters(epsilon=1.5, sigma=2.0, cutoff=3.0, shift=true)
            @test typeof(lj.epsilon) === typeof(1.0u"eV")
            @test typeof(lj.sigma) === typeof(1.0u"Å")
            @test typeof(lj.cutoff) === Float64
            @test typeof(lj.shift) === typeof(1.0u"eV")
        end
    end


    @testset "lj_energy functions tests" begin
        @testset "Basic tests" begin
            # Test at characteristic distances
            epsilon = 1.0u"eV"
            sigma = 1.0u"Å"
            
            # At potential crosses zero
            @test lj_energy(epsilon, sigma, sigma) ≈ 0.0u"eV" atol=1e-10u"eV"
            
            # At the potential minimum
            r_min = 2.0^(1/6) * sigma
            e_min = lj_energy(epsilon, sigma, r_min)
            @test e_min < 0.0u"eV"
            @test e_min < lj_energy(epsilon, sigma, 0.99 * r_min)
            @test e_min < lj_energy(epsilon, sigma, 1.01 * r_min)
            
            # Test repulsive region
            @test lj_energy(epsilon, sigma, 0.8sigma) > 0.0u"eV"
            @test lj_energy(epsilon, sigma, 0.5sigma) > lj_energy(epsilon, sigma, 0.8sigma)
            
            # Test attractive region
            @test lj_energy(epsilon, sigma, 1.5sigma) < 0.0u"eV"
            @test abs(lj_energy(epsilon, sigma, 5.0sigma)) < abs(lj_energy(epsilon, sigma, 2.0sigma))
        end

        @testset "Energy scaling tests" begin
            r = 1.5u"Å"
            sigma = 1.0u"Å"
            
            # Test epsilon scaling
            e1 = lj_energy(1.0u"eV", sigma, r)
            e2 = lj_energy(2.0u"eV", sigma, r)
            @test e2 ≈ 2 * e1
            
            # Test sigma scaling
            r = 3.0u"Å"
            e1 = lj_energy(1.0u"eV", 1.0u"Å", r)
            e2 = lj_energy(1.0u"eV", 2.0u"Å", r)
            @test abs(e2) > abs(e1)
        end

        @testset "lj_energy with LJParameters tests" begin
            # Test basic functionality with default parameters
            lj = LJParameters()
            r = 1.5u"Å"
            e = lj_energy(r, lj)
            @test typeof(e) == typeof(1.0u"eV")
            @test e < 0.0u"eV"
            
            # Test cutoff behavior
            lj_cut = LJParameters(cutoff=2.0, shift=false)
            @test lj_energy(1.9u"Å", lj_cut) != 0.0u"eV"
            @test lj_energy(2.1u"Å", lj_cut) == 0.0u"eV"
            
            # Test shift behavior
            lj_shifted = LJParameters(cutoff=2.0, shift=true)
            @test lj_energy(2.0u"Å", lj_shifted) ≈ 0.0u"eV" atol=1e-10u"eV"
            @test lj_energy(1.5u"Å", lj_shifted) != lj_energy(1.5u"Å", lj_cut)
        end

        @testset "Physical consistency tests" begin
            # Test long-range behavior
            lj = LJParameters(cutoff=Inf)
            @test abs(lj_energy(100.0u"Å", lj)) < 1e-10u"eV"
            
            # Test energy profile shape
            r_values = [0.8, 0.9, 1.0, 1.122, 1.5, 2.0]u"Å"
            energies = [lj_energy(r, lj) for r in r_values]
            
            @test all(diff(energies[1:4]) .< 0u"eV")
            @test all(diff(energies[4:end]) .> 0u"eV")
        end
    
        @testset "Edge cases tests" begin
            # Test with zero epsilon
            lj = LJParameters(epsilon=0.0)
            @test lj_energy(1.5u"Å", lj) == 0.0u"eV"
            
            # Test with very small distances
            lj = LJParameters()
            @test lj_energy(1e-10u"Å", lj) > 1e10u"eV"
            
            # Test with very large distances
            @test abs(lj_energy(1e10u"Å", lj)) < 1e-10u"eV"
            
            # Test with equal sigma and r
            @test lj_energy(1.0u"eV", 2.0u"Å", 2.0u"Å") ≈ 0.0u"eV" atol=1e-10u"eV"
        end
    
        @testset "Numerical precision tests" begin
            epsilon = 1.0u"eV"
            sigma = 1.0u"Å"
            r = 1.122462048309373u"Å"
            
            # Test energy at minimum is close to analytical value
            e_min = lj_energy(epsilon, sigma, r)
            analytical_min = -1.0u"eV"
            @test e_min ≈ analytical_min rtol=1e-3
            
            # Test energy differences at very large distances
            r1 = 10.0u"Å"
            r2 = 10.1u"Å"
            @test abs(lj_energy(epsilon, sigma, r1) - lj_energy(epsilon, sigma, r2)) < 1e-6u"eV"
        end
    
        @testset "Unit consistency tests" begin
            # Test that output maintains correct units
            e = lj_energy(1.0u"eV", 1.0u"Å", 2.0u"Å")
            @test typeof(e) == typeof(1.0u"eV")
        end
    end


    @testset "CompositeLJParameters Tests" begin
        @testset "Direct matrix constructor tests" begin
            # Test 2x2 matrix constructor
            lj_matrix_2x2 = [
                LJParameters(epsilon=1.0) LJParameters(epsilon=2.0);
                LJParameters(epsilon=2.0) LJParameters(epsilon=3.0)
            ]
            comp_lj = CompositeLJParameters{2}(lj_matrix_2x2)
            @test size(comp_lj.lj_param_sets) == (2, 2)
            @test comp_lj.lj_param_sets[1,1].epsilon == 1.0u"eV"
            @test comp_lj.lj_param_sets[2,2].epsilon == 3.0u"eV"
    
            # Test error for wrong matrix size
            @test_throws ArgumentError CompositeLJParameters{3}(lj_matrix_2x2)
        end
    
        @testset "Vector constructor - full matrix tests" begin
            # Test 2x2 case
            ljs_2x2 = [
                LJParameters(epsilon=1.0),  # (1,1)
                LJParameters(epsilon=2.0),  # (1,2)
                LJParameters(epsilon=3.0),  # (2,1)
                LJParameters(epsilon=4.0)   # (2,2)
            ]
            comp_lj = CompositeLJParameters(2, ljs_2x2)
            @test size(comp_lj.lj_param_sets) == (2, 2)
            @test comp_lj.lj_param_sets[1,1].epsilon == 1.0u"eV"
            @test comp_lj.lj_param_sets[2,1].epsilon == 2.0u"eV"
            @test comp_lj.lj_param_sets[1,2].epsilon == 3.0u"eV"
            @test comp_lj.lj_param_sets[2,2].epsilon == 4.0u"eV"
    
            # Test 3x3 case from docstring example
            ljs_3x3 = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]
            comp_lj = CompositeLJParameters(3, ljs_3x3)
            @test size(comp_lj.lj_param_sets) == (3, 3)
            @test comp_lj.lj_param_sets[1,1].epsilon == 11.0u"eV"
            @test comp_lj.lj_param_sets[3,3].epsilon == 33.0u"eV"
            
            # Test error for wrong number of parameters
            @test_throws ArgumentError CompositeLJParameters(2, ljs_3x3)
        end
    
        @testset "Vector constructor - symmetric matrix tests" begin
            # Test 2x2 symmetric case
            ljs_sym_2x2 = [
                LJParameters(epsilon=1.0),  # (1,1)
                LJParameters(epsilon=2.0),  # (1,2) and (2,1)
                LJParameters(epsilon=3.0)   # (2,2)
            ]
            comp_lj = CompositeLJParameters(2, ljs_sym_2x2)
            @test size(comp_lj.lj_param_sets) == (2, 2)
            @test comp_lj.lj_param_sets[1,2] == comp_lj.lj_param_sets[2,1]
            @test comp_lj.lj_param_sets[1,1].epsilon == 1.0u"eV"
            @test comp_lj.lj_param_sets[2,2].epsilon == 3.0u"eV"
    
            # Test 3x3 symmetric case
            ljs_sym_3x3 = [
                LJParameters(epsilon=1.0),  # (1,1)
                LJParameters(epsilon=2.0),  # (1,2) and (2,1)
                LJParameters(epsilon=3.0),  # (1,3) and (3,1)
                LJParameters(epsilon=4.0),  # (2,2)
                LJParameters(epsilon=5.0),  # (2,3) and (3,2)
                LJParameters(epsilon=6.0)   # (3,3)
            ]
            comp_lj = CompositeLJParameters(3, ljs_sym_3x3)
            @test size(comp_lj.lj_param_sets) == (3, 3)
            @test comp_lj.lj_param_sets[1,2] == comp_lj.lj_param_sets[2,1]
            @test comp_lj.lj_param_sets[1,3] == comp_lj.lj_param_sets[3,1]
            @test comp_lj.lj_param_sets[2,3] == comp_lj.lj_param_sets[3,2]
        end
    
        @testset "Equality operator tests" begin
            # Test equal CompositeLJParameters
            ljs1 = [LJParameters(epsilon=e) for e in [1, 2, 3, 4]]
            ljs2 = [LJParameters(epsilon=e) for e in [1, 2, 3, 4]]
            comp_lj1 = CompositeLJParameters(2, ljs1)
            comp_lj2 = CompositeLJParameters(2, ljs2)
            @test comp_lj1 == comp_lj2
    
            # Test unequal CompositeLJParameters
            ljs3 = [LJParameters(epsilon=e) for e in [1, 2, 3, 5]]  # Different last value
            comp_lj3 = CompositeLJParameters(2, ljs3)
            @test comp_lj1 != comp_lj3
        end
    
        @testset "Show method tests" begin
            # Create a test instance
            ljs = [LJParameters(epsilon=e) for e in [1, 2, 2, 3]]
            comp_lj = CompositeLJParameters(2, ljs)
            
            # Capture output
            io = IOBuffer()
            show(io, comp_lj)
            output = String(take!(io))
            
            # Test output format
            @test contains(output, "CompositeLJParameters{2}")
            @test contains(output, "2x2 Matrix{LJParameters}")
            @test contains(output, "lj_param_sets[1, 1]")
            @test contains(output, "LJParameters(1.0 eV")
            
            # Test output for larger matrix
            ljs_3x3 = [LJParameters(epsilon=e) for e in [11, 12, 13, 22, 23, 33]]
            comp_lj_3 = CompositeLJParameters(3, ljs_3x3)
            io = IOBuffer()
            show(io, comp_lj_3)
            output = String(take!(io))
            @test contains(output, "CompositeLJParameters{3}")
            @test contains(output, "3x3 Matrix{LJParameters}")
            @test contains(output, "lj_param_sets[3, 3]")
        end
    
        @testset "Error cases tests" begin
            # Test invalid matrix sizes
            @test_throws ArgumentError CompositeLJParameters(4, [LJParameters() for _ in 1:3])
            @test_throws ArgumentError CompositeLJParameters(3, [LJParameters() for _ in 1:5])
            
            # Test invalid number of parameters for symmetric case
            @test_throws ArgumentError CompositeLJParameters(3, [LJParameters() for _ in 1:4])
            
            # Test attempting to create non-square matrix
            lj_matrix_nonsquare = [
                LJParameters() LJParameters();
                LJParameters() LJParameters();
                LJParameters() LJParameters()
            ]
            @test_throws ArgumentError CompositeLJParameters{2}(lj_matrix_nonsquare)
        end
    end

end