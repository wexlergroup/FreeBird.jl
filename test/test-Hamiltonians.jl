@testset "Hamiltonians Tests" begin

    @testset "Lattice-Gas Hamiltonian" begin
        
        @test typeof(LatticeGasHamiltonian) == DataType

        lattice_gas_hamiltonian = LatticeGasHamiltonian(1u"eV", 2.2u"eV",3.33u"eV")

        @test lattice_gas_hamiltonian |> typeof == LatticeGasHamiltonian

        @test lattice_gas_hamiltonian.adsorption_energy == 1u"eV"
        @test typeof(lattice_gas_hamiltonian.adsorption_energy) == typeof(1.0u"eV")

        @test lattice_gas_hamiltonian.nn_interaction_energy == 2.2u"eV"
        @test typeof(lattice_gas_hamiltonian.nn_interaction_energy) == typeof(1.0u"eV")

        @test lattice_gas_hamiltonian.nnn_interaction_energy == 3.33u"eV"
        @test typeof(lattice_gas_hamiltonian.nnn_interaction_energy) == typeof(1.0u"eV")

    end

    @testset "GenericLatticeHamiltonian" begin

        @test_throws ArgumentError GenericLatticeHamiltonian{2,Float64}(1.0, [1.0])
        @test_throws ArgumentError GenericLatticeHamiltonian{1,Float64}(1.0, [1.0, 2.0])

        h = GenericLatticeHamiltonian(1.0, [2.0, 3.0])
        @test typeof(h) == GenericLatticeHamiltonian{2,Float64}
        @test h.on_site_interaction == 1.0
        @test h.nth_neighbor_interactions[1] == 2.0
        @test h.nth_neighbor_interactions[2] == 3.0

        h = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        @test unit(h.on_site_interaction) == u"eV"
        @test h.on_site_interaction == -0.04u"eV"
        @test all(unit.(h.nth_neighbor_interactions) .== u"eV")
        @test h.nth_neighbor_interactions[1] == -0.01u"eV"
        @test h.nth_neighbor_interactions[2] == -0.0025u"eV"

        io = IOBuffer()
        show(io, h)
        output = String(take!(io))
        @test contains(output, "GenericLatticeHamiltonian")
        @test contains(output, "on_site_interaction")
        @test contains(output, "nth_neighbor_interactions")

    end

    @testset "GenericLatticeHamiltonian" begin

        base_ham = GenericLatticeHamiltonian(-0.04, [-0.01, -0.0025], u"eV")
        
        hams = [deepcopy(base_ham) for _ in 1:4]
        mlham = Hamiltonians.MLatticeHamiltonian(2, hams)
        io = IOBuffer()
        show(io, mlham)
        output = String(take!(io))
        @test size(mlham.Hamiltonians) == (2, 2)
        @test mlham.Hamiltonians[1,1] == base_ham
        @test mlham.Hamiltonians[2,2] == base_ham
        @test contains(output, "MLatticeHamiltonian")
        @test contains(output, "Hamiltonians[1, 1]")
        @test contains(output, "Hamiltonians[2, 2]")

        hams = [deepcopy(base_ham) for _ in 1:3]  # 2x2 symmetric needs 3 elements
        mlham = MLatticeHamiltonian(2, hams)
        @test size(mlham.Hamiltonians) == (2, 2)
        @test mlham.Hamiltonians[1,2] == mlham.Hamiltonians[2,1]

        hams = [deepcopy(base_ham) for _ in 1:2]
        @test_throws ArgumentError MLatticeHamiltonian(2, hams)  # Wrong number of elements

        wrong_matrix = reshape([deepcopy(base_ham) for _ in 1:6], 2, 3)
        @test_throws ArgumentError MLatticeHamiltonian{2,2,typeof(1.0u"eV")}(wrong_matrix)


    end

end 