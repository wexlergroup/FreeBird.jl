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