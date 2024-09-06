@testset "Lattice-Gas Hamiltonian" begin
    @test typeof(LGHamiltonian) == DataType

    lg_h = LGHamiltonian(1u"eV", 2.2u"eV",3.33u"eV")

    @test lg_h |> typeof == LGHamiltonian
    @test lg_h.adsorption_energy == 1u"eV"
    @test typeof(lg_h.adsorption_energy) == typeof(1.0u"eV")
    @test lg_h.nn_interaction_energy == 2.2u"eV"
    @test typeof(lg_h.nn_interaction_energy) == typeof(1.0u"eV")
    @test lg_h.nnn_interaction_energy == 3.33u"eV"
    @test typeof(lg_h.nnn_interaction_energy) == typeof(1.0u"eV")
end 