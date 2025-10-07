@testset "AbstractHamiltonians Tests" begin

    @testset "GenericLatticeHamiltonian property tests" begin
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
            hams_2x2_sym = [GenericLatticeHamiltonian(-0.04, [-0.01*i, -0.0025*i], u"eV") for i in 1:3]
            mlham = MLatticeHamiltonian(2, hams_2x2_sym)
            @test size(mlham.Hamiltonians) == (2, 2)
            # Test symmetry
            @test mlham.Hamiltonians[1,2] === mlham.Hamiltonians[2,1]
    
            # Test 3x3 symmetric matrix (6 elements)
            hams_3x3_sym = [GenericLatticeHamiltonian(-0.04, [-0.01*i, -0.0025*i], u"eV") for i in 1:6]
            mlham = MLatticeHamiltonian(3, hams_3x3_sym)
            @test size(mlham.Hamiltonians) == (3, 3)
            
            # Test symmetry
            @test mlham.Hamiltonians[1,2] === mlham.Hamiltonians[2,1]
            @test mlham.Hamiltonians[1,3] === mlham.Hamiltonians[3,1]
            @test mlham.Hamiltonians[2,3] === mlham.Hamiltonians[3,2]
        end

    end
end