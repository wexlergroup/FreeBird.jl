@testset "AbstractHamiltonians Tests" begin

    @testset "GenericLatticeHamiltonian struct and functions tests" begin

        @testset "Constructor Tests" begin
            # Test error handling for mismatched lengths
            @test_throws ArgumentError GenericLatticeHamiltonian{3,typeof(1.0u"eV")}(1.0u"eV", [1.0, 2.0] .* u"eV")
        end

        @testset "Property Tests" begin
            # Create a test Hamiltonian
            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            
            # Test on-site interaction
            @test ham.on_site_interaction == -0.04u"eV"
            
            # Test nth neighbor interactions
            @test length(ham.nth_neighbor_interactions) == 2
            @test ham.nth_neighbor_interactions[1] == -0.01u"eV"
            @test ham.nth_neighbor_interactions[2] == -0.0025u"eV"
            
            # Test that nth_neighbor_interactions is a StaticVector
            @test ham.nth_neighbor_interactions isa StaticVector
        end
    
        @testset "Unit Consistency" begin
            # Test with different energy units
            ham_ev = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            ham_mev = GenericLatticeHamiltonian(-40.0, [-10.0, -2.5], u"meV")
            
            # Test unit conversion
            @test ham_ev.on_site_interaction ≈ ham_mev.on_site_interaction
            @test ham_ev.nth_neighbor_interactions[1] ≈ ham_mev.nth_neighbor_interactions[1]
            @test ham_ev.nth_neighbor_interactions[2] ≈ ham_mev.nth_neighbor_interactions[2]
        end
    
        @testset "Show Method" begin
            ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
            
            # Test output format
            buf = IOBuffer()
            show(buf, ham)
            output = String(take!(buf))
            
            @test contains(output, "GenericLatticeHamiltonian")
            @test contains(output, "on_site_interaction")
            @test contains(output, "nth_neighbor_interactions")
        end

    end

    @testset "MLatticeHamiltonian struct and functions tests" begin

        @testset "Constructor Tests - Full Matrix" begin
            # Test 2x2 matrix (4 elements)
            hams_2x2 = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
            mlham = MLatticeHamiltonian(2, hams_2x2)
            @test size(mlham.Hamiltonians) == (2, 2)
    
            # Test 3x3 matrix (9 elements)
            hams_3x3 = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:9]
            mlham = MLatticeHamiltonian(3, hams_3x3)
            @test size(mlham.Hamiltonians) == (3, 3)
        end
    
        @testset "Constructor Tests - Symmetric Matrix" begin
            # Test 2x2 symmetric matrix (3 elements)
            hams_2x2_sym = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:3]
            mlham = MLatticeHamiltonian(2, hams_2x2_sym)
            @test size(mlham.Hamiltonians) == (2, 2)
            # Test symmetry
            @test mlham.Hamiltonians[1,2] === mlham.Hamiltonians[2,1]
    
            # Test 3x3 symmetric matrix (6 elements)
            hams_3x3_sym = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:6]
            mlham = MLatticeHamiltonian(3, hams_3x3_sym)
            @test size(mlham.Hamiltonians) == (3, 3)
            
            # Test symmetry
            @test mlham.Hamiltonians[1,2] === mlham.Hamiltonians[2,1]
            @test mlham.Hamiltonians[1,3] === mlham.Hamiltonians[3,1]
            @test mlham.Hamiltonians[2,3] === mlham.Hamiltonians[3,2]
        end

        @testset "Unit Consistency" begin
            # Create Hamiltonians with different units
            hams_ev = [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4]
            hams_mev = [GenericLatticeHamiltonian(-40.0, [-10.0, -2.5], u"meV") for i in 1:4]
            
            mlham_ev = MLatticeHamiltonian(2, hams_ev)
            mlham_mev = MLatticeHamiltonian(2, hams_mev)
    
            # Test unit conversion
            for i in 1:2, j in 1:2
                @test mlham_ev.Hamiltonians[i,j].on_site_interaction ≈ 
                      mlham_mev.Hamiltonians[i,j].on_site_interaction
                @test all(mlham_ev.Hamiltonians[i,j].nth_neighbor_interactions .≈ 
                         mlham_mev.Hamiltonians[i,j].nth_neighbor_interactions)
            end
        end
    
        @testset "Show Method" begin
            mlham = MLatticeHamiltonian(2, [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:4])
            
            # Test output format
            buf = IOBuffer()
            show(buf, mlham)
            output = String(take!(buf))
            
            @test contains(output, "MLatticeHamiltonian")
            @test contains(output, "Hamiltonians[1, 1]")
            @test contains(output, "Hamiltonians[2, 2]")
        end

        @testset "Edge Cases" begin
            # Test error cases
            @test_throws ArgumentError MLatticeHamiltonian{3,2,typeof(1.0u"eV")}([GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:2, j in 1:2])
            @test_throws ArgumentError MLatticeHamiltonian(2, [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:5])
            @test_throws ArgumentError MLatticeHamiltonian(3, [GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV") for i in 1:7])
        end

    end
end